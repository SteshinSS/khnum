#include "simulator/simulator.h"
#include "utilities/debug_utills/debug_prints.h"

#include <iostream>


namespace khnum {
Simulator::Simulator(const std::vector<EmuNetwork> &networks, const std::vector<EmuAndMid> &input_mids,
                     const std::vector<Emu> &measured_isotopes) :
    networks_{networks},
    input_mids_{input_mids},
    measured_isotopes_{measured_isotopes} {}


std::vector<EmuAndMid> Simulator::CalculateMids(const std::vector<Flux>& fluxes) {
    all_known_emus_ = input_mids_;
    fluxes_ = fluxes;
    for (const EmuNetwork &network : networks_) {
        SolveOneNetwork(network);
    }

    return SelectMeasuredMID();
}


// Return vector of simulated MIDs of measured_isotopes
std::vector<EmuAndMid> Simulator::SelectMeasuredMID() {
    std::vector<EmuAndMid> measured_mids;

    for (const Emu &measured_isotope : measured_isotopes_) {
        auto position = find_if(all_known_emus_.begin(),
                                all_known_emus_.end(),
                                [&measured_isotope](const EmuAndMid &emu) {
                                    return emu.emu == measured_isotope;
                                });

        if (position != all_known_emus_.end()) {
            measured_mids.push_back(*position);
        } else {
            throw std::runtime_error("There is a measured isotope which has not computed through metabolic network");
        }
    }

    return measured_mids;
}


void Simulator::SolveOneNetwork(const EmuNetwork &network) {
    const int current_size = FindNetworkSize(network);

    // Solve AX = BY equation
    // See Antoniewitcz 2007

    // EMUs which MIDs are unknown
    // for the X matrix
    std::vector<Emu> unknown_emus;

    // EMUs which MIDs are known
    // for the Y matrix
    std::vector<EmuAndMid> known_emus;

    // So known_emus contains EMUs with known MIDs for this network
    // Whereas the all_known_emus_ has EMUs from other networks

    FillEmuLists(unknown_emus, known_emus, network);

    Matrix A = Matrix::Zero(unknown_emus.size(), unknown_emus.size());
    Matrix B = Matrix::Zero(unknown_emus.size(), known_emus.size());


    Matrix Y = FormYMatrix(known_emus, current_size);
    FillABMatrices(A, B, network, known_emus, unknown_emus);


    Matrix BY = B * Y;
    Matrix X = A.colPivHouseholderQr().solve(BY);

    AppendNewMids(X, unknown_emus, current_size);
}


int Simulator::FindNetworkSize(const EmuNetwork &network) {
    int current_size = 0;
    for (const bool state : network[0].right.emu.atom_states) {
        current_size += static_cast<int>(state);
    }
    return current_size;
}


void Simulator::FillEmuLists(std::vector<Emu> &unknown_emus,
                             std::vector<EmuAndMid> &known_emus,
                             const EmuNetwork &network) {
    // Fills known_emus and unknown_emus
    for (const EmuReaction &reaction : network) {
        if (reaction.left.size() == 1) {
            CheckIsEmuKnown(reaction.left[0].emu, all_known_emus_, known_emus, unknown_emus);
        } else {
            EmuAndMid convolution = ConvolveEmu(reaction.left);
            known_emus.push_back(convolution);
        }

        unknown_emus.push_back(reaction.right.emu);
    }

    // delete repeated emus
    std::sort(known_emus.begin(), known_emus.end());
    known_emus.erase(std::unique(known_emus.begin(), known_emus.end()), known_emus.end());

    std::sort(unknown_emus.begin(), unknown_emus.end());
    unknown_emus.erase(std::unique(unknown_emus.begin(), unknown_emus.end()), unknown_emus.end());
}


void Simulator::CheckIsEmuKnown(const Emu &emu, const std::vector<EmuAndMid> &where_find_emu_list,
                                std::vector<EmuAndMid> &known_emus, std::vector<Emu> &unknown_emus) {
    const Mid *mid = FindMid(emu, where_find_emu_list);
    if (mid) {
        EmuAndMid new_known;
        new_known.emu = emu;
        new_known.mid = *mid;
        known_emus.push_back(new_known);
    } else {
        unknown_emus.push_back(emu);
    }
}


const Mid *Simulator::FindMid(const Emu &emu,
                              const std::vector<EmuAndMid> &mids) {
    auto position = find_if(mids.begin(),
                            mids.end(),
                            [&emu](const EmuAndMid &known_mid) {
                                return known_mid.emu == emu;
                            });
                                                                                              
    if (position == mids.end()) {                                                             
        return nullptr;                                                                       
    } else {                                                                                  
        return &(position->mid);                                                              
    }                                                                                         
}                                                                                             
                                                                                              
                                                                                              
EmuAndMid Simulator::ConvolveEmu(const EmuReactionSide &convolve_reaction) {                  
    // ToDo test function for reaction with 3 and more substrate                              
    EmuAndMid convolution;                                                                    
    convolution.mid = Mid(1, 1.0); // MID = [1.0]                                             
    for (const EmuSubstrate &emu_substrate : convolve_reaction) {                             
        const Emu emu = emu_substrate.emu;                                                    
        convolution.emu.name += emu.name;                                                     
        for (const bool &state : emu.atom_states) {                                           
            convolution.emu.atom_states.push_back(state);                                     
        }
        Mid new_mid = *FindMid(emu, all_known_emus_);
        convolution.mid = convolution.mid * new_mid;
    }

    return convolution;
}


Matrix Simulator::FormYMatrix(const std::vector<EmuAndMid> &known_emus,
                              const int current_size) {
    Matrix Y(known_emus.size(), current_size + 1);
    for (int known_emu_index = 0; known_emu_index < known_emus.size(); ++known_emu_index) {
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            Y(known_emu_index, mass_shift) = known_emus[known_emu_index].mid[mass_shift];
        }
    }
    return Y;
}


void Simulator::FillABMatrices(Matrix &A, Matrix &B,
                               const EmuNetwork &network,
                               const std::vector<EmuAndMid> &known_emus,
                               const std::vector<Emu> &unknown_emus) {
    for (const EmuReaction &reaction : network) {
        EmuSubstrate substrate;
        if (reaction.left.size() > 1) {
            EmuAndMid convolution = ConvolveEmu(reaction.left);
            substrate.emu = convolution.emu;
            substrate.coefficient = 1.0;
        } else {
            substrate = reaction.left[0];
        }


        int position_of_product = FindUnknownEmuPosition(reaction.right.emu, unknown_emus);
        A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes_.at(reaction.id);

        // Return nullptr if there is no substrate.emu in known_emus
        bool is_emu_known = static_cast<bool>(FindMid(substrate.emu, known_emus));
        if (!is_emu_known) {
            int position_of_substrate = FindUnknownEmuPosition(substrate.emu, unknown_emus);
            A(position_of_product, position_of_substrate) += reaction.right.coefficient * fluxes_.at(reaction.id);
        } else {
            int position_of_substrate = FindKnownEmuPosition(substrate.emu, known_emus);
            B(position_of_product, position_of_substrate) += (-substrate.coefficient) * fluxes_.at(reaction.id);
        }

    }
}


int Simulator::FindUnknownEmuPosition(const Emu &emu,
                                      const std::vector<Emu>& unknown_emus) {
    auto position = find(unknown_emus.begin(),
                         unknown_emus.end(),
                         emu);

    return position - unknown_emus.begin();
}


int Simulator::FindKnownEmuPosition(const Emu &emu,
                                    const std::vector<EmuAndMid>& known_emus) {
    auto position = find_if(known_emus.begin(),
                            known_emus.end(),
                            [&emu](const EmuAndMid &known_mid) {
                                return known_mid.emu == emu;
                            });

    return position - known_emus.begin();
}


void Simulator::AppendNewMids(const Matrix &X,
                              const std::vector<Emu> &unknown_emus,
                              const int current_size) {
    for (int previously_unknown_index = 0; previously_unknown_index < unknown_emus.size(); ++previously_unknown_index) {
        EmuAndMid new_known_emu;
        new_known_emu.emu = unknown_emus[previously_unknown_index];
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            new_known_emu.mid.push_back(X(previously_unknown_index, mass_shift));
        }

        all_known_emus_.push_back(new_known_emu);
    }
}
} // namespace khnum
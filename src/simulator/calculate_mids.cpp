#include "calculate_mids.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"


std::vector<EmuAndMid> CalculateMids(const std::vector<Flux>  &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EmuAndMid> all_known_emus,
                                     const std::vector<Emu> &measured_isotopes) {
    for (const EMUNetwork &network : networks) {
        SolveOneNetwork(fluxes, network, all_known_emus);
    }

    return SelectMeasuredMID(all_known_emus, measured_isotopes);
}


int FindNetworkSize(const EMUNetwork &network) {
    int current_size = 0;
    for (const bool state : network[0].right.emu.atom_states) {
        current_size += static_cast<int>(state);
    }
    return current_size;
}


// Return vector of simulated MIDs of measured_isotopes
std::vector<EmuAndMid> SelectMeasuredMID(const std::vector<EmuAndMid> &all_known_emus,
                                         const std::vector<Emu> &measured_isotopes) {
    std::vector<EmuAndMid> measured_mids;

    for (const Emu &measured_isotope : measured_isotopes) {
        auto position = find_if(all_known_emus.begin(),
                                all_known_emus.end(),
                                [&measured_isotope](const EmuAndMid &emu) {
                                    return emu.emu == measured_isotope;
                                });

        if (position != all_known_emus.end()) {
            measured_mids.push_back(*position);
        } else {
            throw std::runtime_error("There is a measured isotope which has not computed through metabolic network");
        }
    }

    return measured_mids;
}


void SolveOneNetwork(const std::vector<Flux> &fluxes,
                     const EMUNetwork &network,
                     std::vector<EmuAndMid> &all_known_emus) {

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
    // Whereas all_known_emus has EMUs from other networks

    FillEMULists(unknown_emus, known_emus, network, all_known_emus);

    Matrix A = Matrix::Zero(unknown_emus.size(), unknown_emus.size());
    Matrix B = Matrix::Zero(unknown_emus.size(), known_emus.size());

    Matrix Y = FormYMatrix(known_emus, current_size);
    FillABMatrices(A, B, network, known_emus, unknown_emus, fluxes, all_known_emus);

    Matrix BY = B * Y;
    Matrix X = A.colPivHouseholderQr().solve(BY);

    AppendNewMids(X,unknown_emus, all_known_emus, current_size);
}

void FillEMULists(std::vector<Emu> &unknown_emus,
                  std::vector<EmuAndMid> &known_emus,
                  const EMUNetwork &network,
                  const std::vector<EmuAndMid> &all_known_emus) {
    // Fills known_emus and unknown_emus
    for (const EmuReaction &reaction : network) {
        if (reaction.left.size() == 1) {
            CheckIsEmuKnown(reaction.left[0].emu, all_known_emus, known_emus, unknown_emus);
        } else {
            EmuAndMid convolution = ConvolveEmu(reaction.left, all_known_emus);
            known_emus.push_back(convolution);
        }
        // NB: we looking for right emu in known_emus, whereas we looking for left emus in the all_known_emus
        // ToDo recall why so
        CheckIsEmuKnown(reaction.right.emu, known_emus, known_emus, unknown_emus);
    }

    // delete repeated emus
    std::sort(known_emus.begin(), known_emus.end());
    known_emus.erase(std::unique(known_emus.begin(), known_emus.end()), known_emus.end());

    std::sort(unknown_emus.begin(), unknown_emus.end());
    unknown_emus.erase(std::unique(unknown_emus.begin(), unknown_emus.end()), unknown_emus.end());
}


const Mid *FindMid(const Emu &emu,
                   const std::vector<EmuAndMid> &known_mids) {
    auto position = find_if(known_mids.begin(),
                            known_mids.end(),
                            [&emu](const EmuAndMid &known_mid) {
                                return known_mid.emu == emu;
                            });

    if (position == known_mids.end()) {
        return nullptr;
    } else {
        return &(position->mid);
    }
}


EmuAndMid ConvolveEmu(const EmuReactionSide &convolve_reaction,
                      const std::vector<EmuAndMid> &known_mids) {
    EmuAndMid convolution;
    convolution.mid = Mid(1, 1.0); // MID = [1.0]
    for (const EmuSubstrate &emu_substrate : convolve_reaction) {
        const Emu emu = emu_substrate.emu;
        convolution.emu.name += emu.name;
        for (const bool &state : emu.atom_states) {
            convolution.emu.atom_states.push_back(state);
        }
        Mid new_mid = *FindMid(emu, known_mids);
        convolution.mid = convolution.mid * new_mid;
    }

    return convolution;
}



Matrix FormYMatrix(const std::vector<EmuAndMid> &known_emus,
                 const int current_size) {
    Matrix Y(known_emus.size(), current_size + 1);
    for (int known_emu_index = 0; known_emu_index < known_emus.size(); ++known_emu_index) {
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            Y(known_emu_index, mass_shift) = known_emus[known_emu_index].mid[mass_shift];
        }
    }
    return Y;
}

void FillABMatrices(Matrix &A, Matrix &B,
                    const EMUNetwork &network,
                    const std::vector<EmuAndMid> &known_emus,
                    const std::vector<Emu> &unknown_emus,
                    const std::vector<Flux> &fluxes,
                    const std::vector<EmuAndMid> &known_mids) {
    for (const EmuReaction &reaction : network) {
        EmuSubstrate substrate;
        if (reaction.left.size() > 1) {
            EmuAndMid convolution = ConvolveEmu(reaction.left, known_mids);
            substrate.emu = convolution.emu;
            substrate.coefficient = 1.0;
        } else {
            substrate = reaction.left[0];
        }

        // return nullptr if there is no substrate.emu in known_emus
        bool is_emu_known = FindMid(substrate.emu, known_emus);
        if (!is_emu_known) {
            // both substrate and product are unknown
            int position_of_substrate = FindUnknownEmuPosition(substrate.emu, unknown_emus);
            int position_of_product = FindUnknownEmuPosition(reaction.right.emu, unknown_emus);
            A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.id);
            A(position_of_product, position_of_substrate) += reaction.right.coefficient * fluxes.at(reaction.id);

        } else {
            // Product is unknown, Substrate is known
            int position_of_substrate = FindKnownEmuPosition(substrate.emu, known_emus);
            int position_of_product = FindUnknownEmuPosition(reaction.right.emu, unknown_emus);

            A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.id);

            B(position_of_product, position_of_substrate) +=
                    (-substrate.coefficient) * fluxes.at(reaction.id);
        }
    }
}



int FindUnknownEmuPosition(const Emu &emu,
                           const std::vector<Emu> unknown_emus) {
    auto position = find(unknown_emus.begin(),
                         unknown_emus.end(),
                         emu);

    return position - unknown_emus.begin();
}

int FindKnownEmuPosition(const Emu &emu,
                         const std::vector<EmuAndMid> known_emus) {
    auto position = find_if(known_emus.begin(),
                            known_emus.end(),
                            [&emu](const EmuAndMid &known_mid) {
                                return known_mid.emu == emu;
                            });

    return position - known_emus.begin();
}


void AppendNewMids(const Matrix &X,
                   const std::vector<Emu> &unknown_emus,
                   std::vector<EmuAndMid> &all_known_emus,
                   const int current_size) {
    for (int previously_unknown_index = 0; previously_unknown_index < unknown_emus.size(); ++previously_unknown_index) {
        EmuAndMid new_known_emu;
        new_known_emu.emu = unknown_emus[previously_unknown_index];
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            new_known_emu.mid.push_back(X(previously_unknown_index, mass_shift));
        }

        all_known_emus.push_back(new_known_emu);
    }
}


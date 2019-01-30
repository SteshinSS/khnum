#include "calculate_mids.h"

#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <exception>


std::vector<EMUandMID> CalculateMids(const std::vector<Flux>  &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EMUandMID> known_mids,
                                     const std::vector<EMU> &measured_isotopes) {
    for (const EMUNetwork &network : networks) {
        SolveOneNetwork(fluxes, network, known_mids);
    }

    return SelectMeasuredMID(known_mids, measured_isotopes);
}


int FindNetworkSize(const EMUNetwork &network) {
    int current_size = 0;
    for (const bool state : network[0].right.emu.atom_states) {
        current_size += static_cast<int>(state);
    }
    return current_size;
}


// Return vector of simulated MIDs of measured_isotopes
std::vector<EMUandMID> SelectMeasuredMID(const std::vector<EMUandMID> &known_mids,
                                         const std::vector<EMU> &measured_isotopes) {
    std::vector<EMUandMID> measured_mids;

    for (const EMU &measured_isotope : measured_isotopes) {
        auto position = find_if(known_mids.begin(),
                                known_mids.end(),
                                [&measured_isotope](const EMUandMID &emu) {
                                    return emu.emu == measured_isotope;
                                });

        if (position != known_mids.end()) {
            measured_mids.push_back(*position);
        } else {
            throw std::runtime_error("There is measured isotope which haven't computed through metabolic network");
        }
    }

    return measured_mids;
}


void SolveOneNetwork(const std::vector<Flux> &fluxes,
                     const EMUNetwork &network,
                     std::vector<EMUandMID> &known_mids) {

    const int current_size = FindNetworkSize(network);

    // Solve AX = BY equation
    // See Antoniewitcz 2007

    // EMUs which MIDs are unknown
    // for the X matrix
    std::vector<EMU> unknown_emus;

    // EMUs which MIDs are known
    // for the Y matrix
    std::vector<EMUandMID> known_emus;

    FillEMULists(unknown_emus, known_emus, network, known_mids);

    Matrix A = Matrix::Zero(unknown_emus.size(), unknown_emus.size());
    Matrix B = Matrix::Zero(unknown_emus.size(), known_emus.size());
    Matrix Y(known_emus.size(), current_size + 1);

    FormYMatrix(Y, known_emus, current_size);
    FormABMatrices(A, B, network, known_emus, unknown_emus, fluxes, known_mids);

    Matrix BY = B * Y;
    Matrix X = A.colPivHouseholderQr().solve(BY);

    AppendNewMIDs(X,unknown_emus, known_mids, current_size);
    return;
}

void FillEMULists(std::vector<EMU> &unknown_emus,
                  std::vector<EMUandMID> &known_emus,
                  const EMUNetwork &network,
                  const std::vector<EMUandMID> &known_mids) {

    // Fills known_emus and unknown_emus
    for (const EMUReaction &reaction : network) {

        // checking the left side
        if (reaction.left.size() == 1) {
            const MID *mid = GetMID(reaction.left[0].emu, known_mids);
            if (mid) {
                EMUandMID new_known;
                new_known.emu = reaction.left[0].emu;
                new_known.mid = *mid;
                known_emus.push_back(new_known);
            } else {
                unknown_emus.push_back(reaction.left[0].emu);
            }
        } else {
            EMUandMID convolution = ConvolveEMU(reaction.left, known_mids);
            known_emus.push_back(convolution);
        }

        // checking the right side
        const MID *mid = GetMID(reaction.right.emu, known_emus);
        if (mid) {
            EMUandMID new_known;
            new_known.emu = reaction.right.emu;
            new_known.mid = *mid;
            known_emus.push_back(new_known);
        } else {
            unknown_emus.push_back(reaction.right.emu);
        }
    }

    // delete repeated emus
    std::sort(known_emus.begin(), known_emus.end());
    known_emus.erase(std::unique(known_emus.begin(), known_emus.end()), known_emus.end());

    std::sort(unknown_emus.begin(), unknown_emus.end());
    unknown_emus.erase(std::unique(unknown_emus.begin(), unknown_emus.end()), unknown_emus.end());

    return;
}


void FormYMatrix(Matrix &Y,
                 const std::vector<EMUandMID> &known_emus,
                 const int current_size) {
    for (int known_emu_index = 0; known_emu_index < known_emus.size(); ++known_emu_index) {
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            Y(known_emu_index, mass_shift) = known_emus[known_emu_index].mid[mass_shift];
        }
    }
}

void FormABMatrices(Matrix &A, Matrix &B,
                    const EMUNetwork &network,
                    const std::vector<EMUandMID> &known_emus,
                    const std::vector<EMU> &unknown_emus,
                    const std::vector<Flux> &fluxes,
                    const std::vector<EMUandMID> &known_mids) {
    for (const EMUReaction &reaction : network) {
        EMUSubstrate substrate;
        if (reaction.left.size() > 1) {
            EMUandMID convolution = ConvolveEMU(reaction.left, known_mids);
            substrate.emu = convolution.emu;
            substrate.coefficient = 1.0;
        } else {
            substrate = reaction.left[0];
        }

        if (!IsEMUKnown(substrate.emu, known_emus)) {
            // they are both unknown
            int position_of_substrate = FindUnknownEMUsPosition(substrate.emu, unknown_emus);
            int position_of_product = FindUnknownEMUsPosition(reaction.right.emu, unknown_emus);
            A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.id);

            // Why does it multiple by product coefficient? Shouldn't it be substrate coefficient?
            A(position_of_product, position_of_substrate) += reaction.right.coefficient * fluxes.at(reaction.id);
        } else {
            // Product is unknown, Substrate is known

            int position_of_substrate = FindKnownEMUsPosition(substrate.emu, known_emus);
            int position_of_product = FindUnknownEMUsPosition(reaction.right.emu, unknown_emus);

            A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.id);

            B(position_of_product, position_of_substrate) +=
                    (-substrate.coefficient) * fluxes.at(reaction.id);
        }

    }
    return;
}


void AppendNewMIDs(const Matrix &X,
                   const std::vector<EMU> &unknown_emus,
                   std::vector<EMUandMID> &known_mids,
                   const int current_size) {

    for (int previously_unknown_index = 0; previously_unknown_index < unknown_emus.size(); ++previously_unknown_index) {
        EMUandMID new_known_emu;
        new_known_emu.emu = unknown_emus[previously_unknown_index];
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            new_known_emu.mid.push_back(X(previously_unknown_index, mass_shift));
        }

        known_mids.push_back(new_known_emu);
    }

}


bool IsEMUKnown(const EMU &emu,
                const std::vector<EMUandMID> known_emus) {
    auto position = find_if(known_emus.begin(),
                            known_emus.end(),
                            [&emu](const EMUandMID &known_mid) {
                                return known_mid.emu == emu;
                            });

    if (position == known_emus.end()) {
        return false;
    } else {
        return true;
    }
}

int FindUnknownEMUsPosition(const EMU &emu,
                            const std::vector<EMU> unknown_emus) {
    auto position = find(unknown_emus.begin(),
                         unknown_emus.end(),
                         emu);

    return position - unknown_emus.begin();
}

int FindKnownEMUsPosition(const EMU &emu,
                          const std::vector<EMUandMID> known_emus) {
    auto position = find_if(known_emus.begin(),
                            known_emus.end(),
                            [&emu](const EMUandMID &known_mid) {
                                return known_mid.emu == emu;
                            });

    return position - known_emus.begin();
}

EMUandMID ConvolveEMU(const EMUReactionSide &convolve_reaction,
                      const std::vector<EMUandMID> &known_mids) {
    EMUandMID convolve_result;
    convolve_result.mid = MID(1, 1.0);
    for (const EMUSubstrate &emu : convolve_reaction) {
        convolve_result.emu.name += emu.emu.name;
        for (const bool &state : emu.emu.atom_states) {
            convolve_result.emu.atom_states.push_back(state);
        }
        MID new_mid = *GetMID(emu.emu, known_mids);
        convolve_result.mid = convolve_result.mid * new_mid;
    }

    return convolve_result;
}

const MID *GetMID(const EMU &emu,
                  const std::vector<EMUandMID> &known_mids) {
    auto position = find_if(known_mids.begin(),
                            known_mids.end(),
                            [&emu](const EMUandMID &known_mid) {
                                return known_mid.emu == emu;
                            });

    if (position == known_mids.end()) {
        return nullptr;
    } else {
        return &(position->mid);
    }
}


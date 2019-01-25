#include "calculate_mids.h"

#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <map>


std::vector<EMUandMID> CalculateMids(const std::map<std::string, Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EMUandMID> known_mids,
                                     const std::vector<EMU> &measured_isotopes) {
    for (const EMUNetwork &network : networks) {
        int current_size = 0;
        for (const bool state : network[0].right.emu.atom_states) {
            current_size += static_cast<int>(state);
        }

        SolveOneNetwork(fluxes, network, known_mids, current_size);
    }

    return known_mids;
}

void SolveOneNetwork(const std::map<std::string, Flux> &fluxes,
                                       const EMUNetwork &network,
                                       std::vector<EMUandMID> &known_mids,
                                       int current_size) {
    // Solve AX = BY equation
    // See Antoniewitcz 2007

    // unknown_emus for X
    std::vector<EMU> unknown_emus;

    // known_emus for Y
    std::vector<EMUandMID> known_emus;

    // creates list of all emus of known and unknown mids
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

    Matrix A = Matrix::Zero(unknown_emus.size(), unknown_emus.size());
    // Matrix X(unknown_emus.size(), current_size + 1);
    Matrix B = Matrix::Zero(unknown_emus.size(), known_emus.size());
    Matrix Y(known_emus.size(), current_size + 1);

    // form Y
    for (int known_emu_index = 0; known_emu_index < known_emus.size(); ++known_emu_index) {
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            Y(known_emu_index, mass_shift) = known_emus[known_emu_index].mid[mass_shift];
        }
    }

    // form A and B

    for (const EMUReaction &reaction : network) {
        if (reaction.left.size() == 1) {
            if (!IsEMUKnown(reaction.left[0].emu, known_emus)) {
                // they are both unknown
                int position_of_substrate = FindUnknownEMUsPosition(reaction.left[0].emu, unknown_emus);
                int position_of_product = FindUnknownEMUsPosition(reaction.right.emu, unknown_emus);
                A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.name);

                // Why does it multiple by product coefficient? Shouldn't it be substrate coefficient?
                A(position_of_product, position_of_substrate) += reaction.right.coefficient * fluxes.at(reaction.name);
            } else {
                // Product is unknown, Substrate is known

                int position_of_substrate = FindKnownEMUsPosition(reaction.left[0].emu, known_emus);
                int position_of_product = FindUnknownEMUsPosition(reaction.right.emu, unknown_emus);

                A(position_of_product, position_of_product) += (-reaction.right.coefficient) * fluxes.at(reaction.name);

                B(position_of_product, position_of_substrate) +=
                        (-reaction.left[0].coefficient) * fluxes.at(reaction.name);
            }
        }
    }

    Matrix BY = B * Y;
    Matrix X = A.colPivHouseholderQr().solve(BY);
    
    for (int previously_unknown_index = 0; previously_unknown_index < unknown_emus.size(); ++previously_unknown_index) {
        EMUandMID new_known_emu;
        new_known_emu.emu = unknown_emus[previously_unknown_index];
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            new_known_emu.mid.push_back(X(previously_unknown_index, mass_shift));
        }

        known_mids.push_back(new_known_emu);
    }

    return;
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


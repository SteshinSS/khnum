#include "calculate_mids.h"

#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"

#include <vector>
#include <algorithm>
#include <iostream>


std::vector<EMUandMID> CalculateMids(const std::vector<Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     const std::vector<EMUandMID> &input_substrates_mids,
                                     const std::vector<EMU> &measured_isotopes) {

}

std::vector<EMUandMID> SolveOneNetwork(const std::vector<Flux> &fluxes,
                                       const EMUNetwork &network,
                                       const std::vector<EMUandMID> &known_mids,
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

    // make unique for known and unknown emus

    Matrix A(unknown_emus.size(), unknown_emus.size());
    // Matrix X(unknown_emus.size(), current_size + 1);
    Matrix B(unknown_emus.size(), known_emus.size());
    Matrix Y(known_emus.size(), current_size + 1);

    // form Y
    for (int known_emu_index = 0; known_emu_index < known_emus.size(); ++known_emu_index) {
        for (int mass_shift = 0; mass_shift < current_size + 1; ++mass_shift) {
            Y(known_emu_index, mass_shift) = known_emus[known_emu_index].mid[mass_shift];
        }
    }

    std::cerr << Y;

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


#include "calculate_mids.h"

#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"

#include <vector>
#include <algorithm>


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
    std::vector<EMU> known_emus;

    for (const EMUReaction &reaction : network) {
        if (reaction.left.size() == 1) {
            if (IsMIDKnown(reaction.left[0].emu, known_mids)) {
                known_emus.push_back(reaction.left[0].emu);
            } else {
                unknown_emus.push_back(reaction.left[0].emu);
            }
        }

    }

    Matrix A(unknown_emus.size(), unknown_emus.size());
    // Matrix X(unknown_emus.size(), current_size + 1);
    Matrix B(unknown_emus.size(), known_emus.size());
    Matrix Y(known_emus.size(), current_size + 1);

    for (int i = 0; i < known_emus.size(); ++i) {

    }

}

bool IsMIDKnown(const EMU &emu,
                const std::vector<EMUandMID> &known_mids) {
    auto position = find_if(known_mids.begin(),
                            known_mids.end(),
                            [&emu](const EMUandMID &known_mid) {
                                return known_mid.emu == emu;
                            });

    if (position == known_mids.end()) {
        return false;
    } else {
        return true;
    }
}


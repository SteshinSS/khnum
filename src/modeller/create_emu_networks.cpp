#include "create_emu_networks.h"

#include "../utilities/EMU.h"
#include "../utilities/MID.h"

#include <vector>
#include <queue>
#include <algorithm>

std::vector<EMUNetwork> CreateEMUNetworks(const std::vector<EMUReaction> &reactions,
                                          const std::vector<EMU> &input_emu_list,
                                          const std::vector<EMU> &measured_isotopes) {
    int max_size = FindTheLargestEMUSize(reactions);

    // emu_networks[i] contains EMU network of i + 1 size
    std::vector<EMUNetwork> emu_networks(max_size);

    // dfs queue
    std::queue<EMU> emus_to_check;
    for (const EMU &measured_isotope : measured_isotopes) {
        emus_to_check.push(measured_isotope);
    }

    // visited
    std::vector<EMU> already_checked_emus;
    for (const EMU &emu : input_emu_list) {
        already_checked_emus.push_back(emu);
    }

    while (!emus_to_check.empty()) {
        EMU next_emu = emus_to_check.front();
        emus_to_check.pop();
        int emu_size = GetEMUSize(next_emu);

        if(IsEMUAlreadyChecked(next_emu, already_checked_emus)) {
            continue;
        }

        for (const EMUReaction &emu_reaction : reactions) {
            if (emu_reaction.right.emu == next_emu) {

                emu_networks[emu_size - 1].push_back(emu_reaction);
                for (const EMUSubstrate &emu_substrate : emu_reaction.left) {
                    if (!IsEMUAlreadyChecked(emu_substrate.emu, already_checked_emus)) {
                        emus_to_check.push(emu_substrate.emu);
                    }
                }
            }
        }
        already_checked_emus.push_back(next_emu);
    }

    return emu_networks;
}

int FindTheLargestEMUSize (const std::vector<EMUReaction> &reactions) {
    int max_size = -1;
    for (const EMUReaction &reaction : reactions) {
        for (const EMUSubstrate &emu_substrate : reaction.left) {
            max_size = std::max(max_size, GetEMUSize(emu_substrate.emu));
        }
        max_size = std::max(max_size, GetEMUSize(reaction.right.emu));
    }
    return max_size;
}

bool IsEMUAlreadyChecked (const EMU &emu, const std::vector<EMU> &already_checked_emus) {
    auto emu_position = find(already_checked_emus.begin(),
                             already_checked_emus.end(),
                             emu);
    return emu_position != already_checked_emus.end();
}

int GetEMUSize(const EMU &emu) {
    int size = 0;
    for (const bool &state : emu.atom_states) {
        if (state) {
            ++size;
        }
    }

    return size;
}
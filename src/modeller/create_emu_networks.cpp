#include "create_emu_networks.h"

#include <vector>
#include <queue>
#include <algorithm>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"


std::vector<EMUNetwork> CreateEMUNetworks(const std::vector<EMUReaction> &reactions,
                                          const std::vector<Emu> &input_emu_list,
                                          const std::vector<Emu> &measured_isotopes) {
    int max_size = FindTheLargestEMUSize(reactions);

    // emu_networks[i] contains Emu network of i + 1 size
    std::vector<EMUNetwork> emu_networks(max_size);

    // dfs queue
    std::queue<Emu> emus_to_check;
    for (const Emu &measured_isotope : measured_isotopes) {
        emus_to_check.push(measured_isotope);
    }

    // visited
    std::vector<Emu> already_checked_emus;
    for (const Emu &emu : input_emu_list) {
        already_checked_emus.push_back(emu);
    }

    while (!emus_to_check.empty()) {
        Emu next_emu = emus_to_check.front();
        emus_to_check.pop();
        int emu_size = GetEMUSize(next_emu);

        if (IsEMUAlreadyChecked(next_emu, already_checked_emus)) {
            continue;
        }

        for (const EMUReaction &emu_reaction : reactions) {
            if (emu_reaction.right.emu == next_emu) {

                emu_networks[emu_size - 1].push_back(emu_reaction);
                for (const EmuSubstrate &emu_substrate : emu_reaction.left) {
                    if (!IsEMUAlreadyChecked(emu_substrate.emu, already_checked_emus)) {
                        emus_to_check.push(emu_substrate.emu);
                    }
                }
            }
        }
        already_checked_emus.push_back(next_emu);
    }

    // remove empty networks

    emu_networks.erase(
            std::remove_if(emu_networks.begin(), emu_networks.end(), [](const EMUNetwork &network) {
                return network.empty();
            }),
            emu_networks.end());

    return emu_networks;
}

int FindTheLargestEMUSize(const std::vector<EMUReaction> &reactions) {
    int max_size = -1;
    for (const EMUReaction &reaction : reactions) {
        for (const EmuSubstrate &emu_substrate : reaction.left) {
            max_size = std::max(max_size, GetEMUSize(emu_substrate.emu));
        }
        max_size = std::max(max_size, GetEMUSize(reaction.right.emu));
    }
    return max_size;
}

bool IsEMUAlreadyChecked(const Emu &emu, const std::vector<Emu> &already_checked_emus) {
    auto emu_position = find(already_checked_emus.begin(),
                             already_checked_emus.end(),
                             emu);
    return emu_position != already_checked_emus.end();
}

int GetEMUSize(const Emu &emu) {
    int size = 0;
    for (const bool &state : emu.atom_states) {
        if (state) {
            ++size;
        }
    }

    return size;
}

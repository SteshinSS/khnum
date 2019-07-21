#include "create_emu_networks.h"

#include <vector>
#include <queue>
#include <algorithm>
#include <iostream>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/debug_utills/debug_prints.h"


namespace khnum {
namespace modelling_utills {
std::vector<EmuNetwork> CreateEmuNetworks(const std::vector<EmuReaction> &reactions,
                                          const std::vector<Emu> &input_emu_list,
                                          const std::vector<Emu> &measured_isotopes) {
    int max_size = FindLargestEmuSize(reactions);

    // emu_networks[i] contains Emu network of i + 1 size
    std::vector<EmuNetwork> emu_networks(max_size);

    // dfs queue
    std::queue<Emu> emus_to_check;
    for (const Emu &measured_isotope : measured_isotopes) {
        emus_to_check.push(measured_isotope);
    }

    std::vector<Emu> already_checked_emus;
    for (const Emu &emu : input_emu_list) {
        already_checked_emus.push_back(emu);
    }

    while (!emus_to_check.empty()) {
        Emu next_emu = emus_to_check.front();
        emus_to_check.pop();
        int emu_size = GetEmuSize(next_emu);

        if (IsEmuAlreadyChecked(next_emu, already_checked_emus)) {
            continue;
        }

        for (const EmuReaction &emu_reaction : reactions) {
            if (emu_reaction.right.emu == next_emu) {
                emu_networks[emu_size - 1].push_back(emu_reaction);
                for (const EmuSubstrate &emu_substrate : emu_reaction.left) {
                    if (!IsEmuAlreadyChecked(emu_substrate.emu, already_checked_emus)) {
                        emus_to_check.push(emu_substrate.emu);
                    }
                }
            }
        }
        already_checked_emus.push_back(next_emu);
    }

    // remove empty networks

    emu_networks.erase(
        std::remove_if(emu_networks.begin(), emu_networks.end(), [](const EmuNetwork &network) {
            return network.empty();
        }),
        emu_networks.end());

    return emu_networks;
}


int FindLargestEmuSize(const std::vector<EmuReaction> &reactions) {
    int max_size = -1;
    for (const EmuReaction &reaction : reactions) {
        for (const EmuSubstrate &emu_substrate : reaction.left) {
            max_size = std::max(max_size, GetEmuSize(emu_substrate.emu));
        }
        max_size = std::max(max_size, GetEmuSize(reaction.right.emu));
    }
    return max_size;
}


int GetEmuSize(const Emu &emu) {
    int size = 0;
    for (const bool &state : emu.atom_states) {
        if (state) {
            ++size;
        }
    }

    return size;
}


bool IsEmuAlreadyChecked(const Emu &emu, const std::vector<Emu> &already_checked_emus) {
    auto emu_position = find(already_checked_emus.begin(),
                             already_checked_emus.end(),
                             emu);
    return emu_position != already_checked_emus.end();
}
} // namespace modelling_utills
} // namespace khnum
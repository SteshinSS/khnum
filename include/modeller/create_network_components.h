#pragma once

#include <vector>

#include "utilities/emu.h"

namespace khnum {
    namespace modelling_utills {
        std::vector<EmuNetwork> CreateNetworkComponents(const EmuNetwork &network);

        void dfs1(int v, const EmuNetwork &network, std::vector<char> &visited, std::vector<size_t> &ordered_reactions, const std::vector<Emu> &all_emus);
        void dfs2(int v, const EmuNetwork &emu_network, std::vector<char> &visited, std::vector<size_t> &component, const std::vector<Emu> &all_emus);
    }
}
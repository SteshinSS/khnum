#pragma once

#include <vector>

#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"


namespace khnum {
struct FluxAndCoefficient {
    int id;
    double coefficient;
};

struct FluxCombination {
    int i, j;
    std::vector<FluxAndCoefficient> fluxes;
};


bool compare(const FluxCombination& lhs, const FluxCombination& rhs);


// Yi contains of emus from the mids_Yi_[network][position]
struct PositionOfKnownEmu {
    int network;
    int position;
};

struct Convolution {
    std::vector<PositionOfKnownEmu> elements;
    int flux_id;
};


struct NetworkEmu {
    Emu emu;
    int network;
    int order_in_X;
    bool is_usefull;
    int order_in_usefull_emus;
};

struct FinalEmu {
    Emu emu;
    int network;
    int order_in_X;
    int position_in_result;
};

bool operator==(const PositionOfKnownEmu& lhs, const PositionOfKnownEmu& rhs);
} // namespace khnum
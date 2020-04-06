#pragma once

#include <vector>

#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"


namespace khnum {
struct FluxAndCoefficient {
    int id;
    double coefficient;
};

struct SparseElement {
    size_t i, j;
    double value;
};

struct FluxCombination {
    size_t i, j;
    std::vector<FluxAndCoefficient> fluxes;
};


bool compare(const FluxCombination& lhs, const FluxCombination& rhs);


// Yi contains of emus from the mids_Yi_[network][position]
struct PositionOfSavedEmu {
    int network;
    int position;
};

struct Convolution {
    std::vector<PositionOfSavedEmu> elements;
    int flux_id;
};

bool operator==(const Convolution& lhs, const Convolution& rhs);


struct NetworkEmu {
    Emu emu;
    int network;
    int order_in_X;
    bool is_usefull;
    int order_in_usefull_emus;
};

struct FinalEmu {
    Emu emu;
    int order_in_X;
    int position_in_result;
    Matrix correction_matrix;
};

bool operator==(const PositionOfSavedEmu& lhs, const PositionOfSavedEmu& rhs);
} // namespace khnum
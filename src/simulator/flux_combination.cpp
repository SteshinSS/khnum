#include "simulator/flux_combination.h"

#include <tuple>


namespace khnum {


bool operator==(const Convolution& lhs, const Convolution& rhs) {
    if (lhs.flux_id == rhs.flux_id) {
        return lhs.elements == rhs.elements;
    } else {
        return false;
    }
}

bool compare(const FluxCombination &lhs, const FluxCombination &rhs) {
    return std::tie(lhs.i, lhs.j) < std::tie(rhs.i, rhs.j);
}

bool operator==(const PositionOfSavedEmu& lhs, const PositionOfSavedEmu& rhs) {
    return std::tie(lhs.network, lhs.position) == std::tie(rhs.network, rhs.position);
}
} // namespace khnum
#include "simulator/flux_combination.h"

#include <tuple>


namespace khnum {
bool compare(const FluxCombination &lhs, const FluxCombination &rhs) {
    return std::tie(lhs.i, lhs.j) < std::tie(rhs.i, rhs.j);
}

bool operator==(const PositionOfKnownEmu& lhs, const PositionOfKnownEmu& rhs) {
    return std::tie(lhs.network, lhs.position) == std::tie(rhs.network, rhs.position);
}
} // namespace khnum
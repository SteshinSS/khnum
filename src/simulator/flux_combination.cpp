#include "simulator/flux_combination.h"

#include <tuple>


namespace khnum {
bool compare(const FluxCombination &lhs, const FluxCombination &rhs) {
    return std::tie(lhs.i, lhs.j) < std::tie(rhs.i, rhs.j);
}

} // namespace khnum
#include "reaction_struct.h"
#include <tuple>

bool operator==(const Substrate &lhs, const Substrate &rhs) {
    return std::tie(lhs.name, lhs.formula, lhs.coefficient) == std::tie(rhs.name, rhs.formula, rhs.coefficient);
}
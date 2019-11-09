#include "utilities/reaction.h"

#include <tuple>


namespace khnum {
bool operator==(const Substrate &lhs, const Substrate &rhs) {
    if (std::tie(lhs.name, lhs.formula) == std::tie(rhs.name, rhs.formula)) {
        double epsilon = 0.00001;
        if ((rhs.substrate_coefficient_ - epsilon < lhs.substrate_coefficient_) &&
            (lhs.substrate_coefficient_ - epsilon < rhs.substrate_coefficient_)) {
            if ((rhs.atom_coefficient_ - epsilon < lhs.atom_coefficient_) &&
                (lhs.atom_coefficient_ - epsilon < rhs.atom_coefficient_)) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}
} // namespace khnum
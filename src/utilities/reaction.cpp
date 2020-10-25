#include "utilities/reaction.h"

#include <tuple>


namespace khnum {
bool operator==(const Substrate &lhs, const Substrate &rhs) {
    if ((lhs.name == rhs.name) && (lhs.id == rhs.id)) {
        double epsilon = 0.00001;
        if ((rhs.substrate_coefficient_ - epsilon < lhs.substrate_coefficient_) &&
            (lhs.substrate_coefficient_ - epsilon < rhs.substrate_coefficient_)) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}
} // namespace khnum
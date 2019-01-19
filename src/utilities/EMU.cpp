#include "EMU.h"
#include <tuple>

bool operator<(const EMU &lhs, const EMU &rhs) {
    return tie(lhs.name, lhs.atom_states) < tie(rhs.name, rhs.atom_states);
}

bool operator==(const EMU &lhs, const EMU &rhs) {
    return !(lhs < rhs) && !(rhs < lhs);
}


bool operator==(EMUSubstrate const &lhs, EMUSubstrate const &rhs) {
    if (lhs.emu == rhs.emu) {
        return lhs.coefficient == rhs.coefficient;
    } else {
        return false;
    }
}

bool operator!=(EMUSubstrate const &lhs, EMUSubstrate const &rhs) {
    return !(lhs == rhs);
}


bool operator==(const EMUReaction &lhs, const EMUReaction &rhs) {
    if (lhs.name != rhs.name) {
        return false;
    } else {
        if (lhs.right != rhs.right) {
            return false;
        } else {
            return lhs.left == rhs.left;
        }
    }
}
#include "utilities/Emu.h"

#include <tuple>

bool operator<(const Emu &lhs, const Emu &rhs) {
    return tie(lhs.name, lhs.atom_states) < tie(rhs.name, rhs.atom_states);
}

bool operator==(const Emu &lhs, const Emu &rhs) {
    return !(lhs < rhs) && !(rhs < lhs);
}


bool operator==(EmuSubstrate const &lhs, EmuSubstrate const &rhs) {
    if (lhs.emu == rhs.emu) {
        return lhs.coefficient == rhs.coefficient;
    } else {
        return false;
    }
}

bool operator!=(EmuSubstrate const &lhs, EmuSubstrate const &rhs) {
    return !(lhs == rhs);
}


bool operator==(const EMUReaction &lhs, const EMUReaction &rhs) {
    if (lhs.id != rhs.id) {
        return false;
    } else {
        if (lhs.right != rhs.right) {
            return false;
        } else {
            return lhs.left == rhs.left;
        }
    }
}
#include "utilities/emu.h"

#include <tuple>


namespace khnum {

// need this for std containers
bool comparator(const Emu& lhs, const Emu& rhs) {
    return tie(lhs.name, lhs.atom_states) < tie(rhs.name, rhs.atom_states);
}

bool operator<(const Emu &lhs, const Emu &rhs) {
    return tie(lhs.name, lhs.atom_states) < tie(rhs.name, rhs.atom_states);
}


bool operator==(const Emu &lhs, const Emu &rhs) {
    return !(lhs < rhs) && !(rhs < lhs);
}


bool operator==(EmuSubstrate const &lhs, EmuSubstrate const &rhs) {
    const double epsilon = 0.0001;
    if (lhs.emu == rhs.emu) {
        bool is_same = ((lhs.coefficient - epsilon) < rhs.coefficient) &&
                       ((rhs.coefficient - epsilon) < lhs.coefficient);
        return is_same;
    } else {
        return false;
    }
}


bool operator!=(EmuSubstrate const &lhs, EmuSubstrate const &rhs) {
    return !(lhs == rhs);
}


bool operator==(const EmuReaction &lhs, const EmuReaction &rhs) {
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
} //namespace khnum
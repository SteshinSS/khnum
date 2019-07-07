#ifndef CFLEX_EMU_H
#define CFLEX_EMU_H

#include <vector>
#include <string>

using AtomStates = std::vector<char>; // like a vector<bool> but iterable
using EMUcoefficient = double;

struct Emu {
    std::string name;
    AtomStates atom_states;
};

// need this for stl containers
bool operator<(const Emu &lhs, const Emu &rhs);
bool operator==(const Emu &lhs, const Emu &rhs);

// uses for Emu reactions
struct EMUSubstrate {
    Emu emu;
    EMUcoefficient coefficient;
};

// need this for stl containers
bool operator==(EMUSubstrate const &lhs, EMUSubstrate const &rhs);
bool operator!=(EMUSubstrate const &lhs, EMUSubstrate const &rhs);

using EMUReactionSide = std::vector<EMUSubstrate>;

struct EMUReaction {
    int id;
    EMUReactionSide left;
    EMUSubstrate right;
};

bool operator==(const EMUReaction &lhs, const EMUReaction &rhs);
using EMUNetwork = std::vector<EMUReaction>;

#endif //CFLEX_EMU_H

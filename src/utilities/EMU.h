#ifndef CFLEX_EMU_H
#define CFLEX_EMU_H

#include <vector>
#include <string>

using AtomStates = std::vector<char>; // there are other available options as vector<char>

using EMUcoefficient = double;

struct EMU {
    std::string name;
    AtomStates atom_states;
};

// need this for stl containers
bool operator<(const EMU &lhs, const EMU &rhs);

bool operator==(const EMU &lhs, const EMU &rhs);

// uses for EMU reactions
struct EMUSubstrate {
    EMU emu;
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

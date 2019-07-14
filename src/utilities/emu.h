#pragma once

#include <vector>
#include <string>

using AtomStates = std::vector<char>; // like a vector<bool> but iterable

struct Emu {
    std::string name;
    AtomStates atom_states;
};


using EmuCoefficient = double;

struct EmuSubstrate {
    Emu emu;
    EmuCoefficient coefficient;
};


using EmuReactionSide = std::vector<EmuSubstrate>;

struct EmuReaction {
    int id;
    EmuReactionSide left;
    EmuSubstrate right;
};

using EmuNetwork = std::vector<EmuReaction>;

// need this for stl containers
bool operator<(const Emu &lhs, const Emu &rhs);

bool operator==(const Emu &lhs, const Emu &rhs);

bool operator==(EmuSubstrate const &lhs, EmuSubstrate const &rhs);

bool operator!=(EmuSubstrate const &lhs, EmuSubstrate const &rhs);

bool operator==(const EmuReaction &lhs, const EmuReaction &rhs);
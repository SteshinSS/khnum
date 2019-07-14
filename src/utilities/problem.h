#ifndef CFLEX_PROBLEM_H
#define CFLEX_PROBLEM_H

#include "math_utilites.h"
#include "input_substrate.h"
#include "MID.h"
#include "Emu.h"
#include "measurements_struct.h"
#include "reaction_struct.h"

struct Problem {
    std::vector<Reaction> reactions;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<EMUNetwork> networks;
    std::vector<EmuAndMid> input_mids;
    std::vector<Measurement> measurements;
    int measurements_count;

};

#endif //CFLEX_PROBLEM_H

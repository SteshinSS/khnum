//
// Created by sema on 28.01.19.
//

#ifndef CFLEX_OBJECTIVE_PARAMETERS_H
#define CFLEX_OBJECTIVE_PARAMETERS_H

#include "math_utilites.h"
#include "input_substrate.h"
#include "MID.h"
#include "Emu.h"
#include "measurements_struct.h"
#include "reaction_struct.h"

struct ObjectiveParameters {
    Matrix *nullspace;
    std::vector<EMUNetwork> *networks;
    std::vector<EmuAndMid> *input_mids;
    std::vector<Emu> *measured_isotopes;
    std::vector<Measurement> *measurements;
    std::vector<Reaction> *reactions;
};

#endif //CFLEX_OBJECTIVE_PARAMETERS_H

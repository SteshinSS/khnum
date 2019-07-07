#ifndef CFLEX_MEASUREMENTS_STRUCT_H
#define CFLEX_MEASUREMENTS_STRUCT_H

#include "Emu.h"
#include "MID.h"

#include <vector>

using Errors = std::vector<double>;

struct Measurement {
    Emu emu;
    MID mid;
    Errors errors;
};

#endif //CFLEX_MEASUREMENTS_STRUCT_H

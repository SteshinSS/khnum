#ifndef CFLEX_MEASUREMENTS_STRUCT_H
#define CFLEX_MEASUREMENTS_STRUCT_H

#include "../utilities/EMU.h"
#include "../utilities/MID.h"

#include <vector>

using Errors = std::vector<double>;

struct Measurement {
    EMU emu;
    MID mid;
    Errors errors;
};

#endif //CFLEX_MEASUREMENTS_STRUCT_H

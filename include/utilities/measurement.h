#pragma once

#include "emu.h"
#include "emu_and_mid.h"
#include "utilities/matrix.h"

#include <vector>


namespace khnum {
using Errors = std::vector<double>;

struct Measurement {
    Emu emu;
    Mid mid;
    Errors errors;
    Matrix correction_matrix;
};
} // namespace khnum
#pragma once

#include "emu.h"
#include "emu_and_mid.h"

#include <vector>

using Errors = std::vector<double>;

struct Measurement {
    Emu emu;
    Mid mid;
    Errors errors;
};
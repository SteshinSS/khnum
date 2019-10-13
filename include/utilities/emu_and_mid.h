#pragma once

#include "utilities/emu.h"

#include <vector>


namespace khnum {
using Mid = std::vector<double>;

struct EmuAndMid {
    Emu emu;
    Mid mid;
};

// convolution
Mid operator*(const Mid &lhs, const Mid &rhs);

// need this for stl containers
bool operator==(const EmuAndMid &lhs, const EmuAndMid &rhs);

bool operator<(const EmuAndMid &lhs, const EmuAndMid &rhs);
} //namespace khnum
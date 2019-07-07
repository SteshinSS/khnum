#ifndef CFLEX_MID_H
#define CFLEX_MID_H

#include "Emu.h"

#include <vector>

using MID = std::vector<double>;

struct EmuAndMid {
  Emu emu;
  MID mid;
};


MID operator*(const MID &lhs, const MID &rhs);
bool operator==(const EmuAndMid &lhs, const EmuAndMid &rhs);

bool operator<(const EmuAndMid &lhs, const EmuAndMid &rhs);
#endif //CFLEX_MID_H

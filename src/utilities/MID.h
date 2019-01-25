#ifndef CFLEX_MID_H
#define CFLEX_MID_H

#include "EMU.h"

#include <vector>

using MID = std::vector<double>;

struct EMUandMID {
  EMU emu;
  MID mid;
};


MID operator*(const MID &lhs, const MID &rhs);
bool operator==(const EMUandMID &lhs, const EMUandMID &rhs);

bool operator<(const EMUandMID &lhs, const EMUandMID &rhs);
#endif //CFLEX_MID_H

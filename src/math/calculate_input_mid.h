#ifndef CFLEX_CALCULATE_INPUT_MID_H
#define CFLEX_CALCULATE_INPUT_MID_H

#include "../utilities/MID.h"
#include "../utilities/input_substrate.h"

#include <vector>

std::vector<EMUandMID> CalculateInputMid(const std::vector<InputSubstrate> &input_substrates,
                                         const std::vector<EMU> &input_emus);


EMUandMID CalculateOneMid(const InputSubstrate &input_substrate,
                          const EMU &input_emu);

#endif //CFLEX_CALCULATE_INPUT_MID_H

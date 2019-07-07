#ifndef CFLEX_CALCULATE_INPUT_MID_H
#define CFLEX_CALCULATE_INPUT_MID_H

#include "MID.h"
#include "input_substrate.h"

#include <vector>

std::vector<EmuAndMid> CalculateInputMid(const std::vector<InputSubstrate> &input_substrates,
                                         const std::vector<Emu> &input_emus);


EmuAndMid CalculateOneMid(const InputSubstrate &input_substrate,
                          const Emu &input_emu);

#endif //CFLEX_CALCULATE_INPUT_MID_H

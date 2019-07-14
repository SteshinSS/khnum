#pragma once

#include <vector>

#include "utilities/emu_and_mid.h"
#include "utilities/input_substrate.h"


std::vector<EmuAndMid> CalculateInputMid(const std::vector<InputSubstrate> &input_substrates,
                                         const std::vector<Emu> &input_emus);


EmuAndMid CalculateOneMid(const InputSubstrate &input_substrate,
                          const Emu &input_emu);
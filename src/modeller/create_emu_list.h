#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/input_substrate.h"


std::vector<Emu> CreateInputEMUList(const std::vector<EMUReaction> &reactions,
                                    const std::vector<InputSubstrate> &input_substrates);


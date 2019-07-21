#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/input_substrate.h"


namespace khnum {
namespace modelling_utills {
std::vector<Emu> CreateInputEmuList(const std::vector<EmuReaction> &reactions,
                                    const std::vector<InputSubstrate> &input_substrates);
} // namespace modelling_utills
} // namespace khnum
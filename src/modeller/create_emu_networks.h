#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"


namespace khnum {
namespace modelling_utills {
std::vector<EmuNetwork> CreateEmuNetworks(const std::vector<EmuReaction> &reactions,
                                          const std::vector<Emu> &input_emu_list,
                                          const std::vector<Emu> &measured_isotopes);

int FindLargestEmuSize(const std::vector<EmuReaction> &reactions);

int GetEmuSize(const Emu &emu);

bool IsEmuAlreadyChecked(const Emu &emu, const std::vector<Emu> &already_checked_emus);
} // namespace modelling_utills
} // namespace khnum
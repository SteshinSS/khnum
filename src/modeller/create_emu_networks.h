#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"


std::vector<EmuNetwork> CreateEMUNetworks(const std::vector<EmuReaction> &reactions,
                                          const std::vector<Emu> &input_emu_list,
                                          const std::vector<Emu> &measured_isotopes);


int FindTheLargestEMUSize(const std::vector<EmuReaction> &reactions);

bool IsEMUAlreadyChecked(const Emu &emu, const std::vector<Emu> &already_checked_emus);

int GetEMUSize(const Emu &emu);

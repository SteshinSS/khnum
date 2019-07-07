#ifndef CFLEX_CREATE_EMU_NETWORKS_H
#define CFLEX_CREATE_EMU_NETWORKS_H

#include "Emu.h"
#include "MID.h"

#include <vector>

std::vector<EMUNetwork> CreateEMUNetworks(const std::vector<EMUReaction> &reactions,
                                          const std::vector<Emu> &input_emu_list,
                                          const std::vector<Emu> &measured_isotopes);


int FindTheLargestEMUSize (const std::vector<EMUReaction> &reactions);
bool IsEMUAlreadyChecked (const Emu &emu, const std::vector<Emu> &already_checked_emus);
int GetEMUSize(const Emu &emu);
#endif //CFLEX_CREATE_EMU_NETWORKS_H

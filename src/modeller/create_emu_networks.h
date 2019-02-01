#ifndef CFLEX_CREATE_EMU_NETWORKS_H
#define CFLEX_CREATE_EMU_NETWORKS_H

#include "EMU.h"
#include "MID.h"

#include <vector>

std::vector<EMUNetwork> CreateEMUNetworks(const std::vector<EMUReaction> &reactions,
                                          const std::vector<EMU> &input_emu_list,
                                          const std::vector<EMU> &measured_isotopes);


int FindTheLargestEMUSize (const std::vector<EMUReaction> &reactions);
bool IsEMUAlreadyChecked (const EMU &emu, const std::vector<EMU> &already_checked_emus);
int GetEMUSize(const EMU &emu);
#endif //CFLEX_CREATE_EMU_NETWORKS_H

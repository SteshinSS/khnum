#ifndef CFLEX_CREATE_EMU_LIST_H
#define CFLEX_CREATE_EMU_LIST_H

#include "Emu.h"
#include "input_substrate.h"

#include <vector>

std::vector<Emu> CreateInputEMUList(const std::vector<EMUReaction> &reactions,
                                    const std::vector<InputSubstrate> &input_substrates);

#endif //CFLEX_CREATE_EMU_LIST_H

//
// Created by Олег on 2019-01-07.
//

#ifndef CFLEX_DEBUG_UTILITES_H
#define CFLEX_DEBUG_UTILITES_H

#include "EMU.h"

const std::string empty_string;

void ShowEMU(const EMU &emu, const std::string &what = empty_string);

void ShowEMUReaction(const EMUReaction &emu_reaction, const std::string &what = empty_string);

#endif //CFLEX_DEBUG_UTILITES_H

#ifndef CFLEX_CREATE_EMU_REACTIONS_H
#define CFLEX_CREATE_EMU_REACTIONS_H

#include "reaction_struct.h"
#include "Emu.h"

#include <vector>
#include <string>
#include <queue>
#include <set>

std::vector<EMUReaction> CreateAllEMUReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<Emu> &observable_emus);


std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const Emu &emu);


std::vector<EMUReaction> CreateNewEMUReactions(const Reaction &reaction,
                                               const Emu &produced_emu);


std::vector<EMUReaction> CreateAllEMUReactions(const Reaction &reaction,
                                               const Emu &produced_emu);

EMUReaction CreateOneEMUReaction(const Reaction &reaction,
                                 const Substrate &produced_emu_substrate,
                                 const Emu &produced_emu);

std::vector<EMUReaction> SelectUniqueEMUReactions(const std::vector<EMUReaction> &all_new_emu_reactions);

void AddNewEMUInQueue(std::queue<Emu> *queue,
                      const std::set<Emu> &emu_ignore_list,
                      const EmuReactionSide &reaction_side);

#endif //CFLEX_CREATE_EMU_REACTIONS_H

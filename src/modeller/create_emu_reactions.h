#ifndef CFLEX_CREATE_EMU_REACTIONS_H
#define CFLEX_CREATE_EMU_REACTIONS_H

#include "../utilities/reaction_struct.h"
#include "../utilities/EMU.h"

#include <vector>
#include <string>
#include <queue>
#include <set>

std::vector<EMUReaction> CreateAllEMUReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<EMU> &observable_emus);


std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const EMU &emu);


std::vector<EMUReaction> CreateNewEMUReactions(const Reaction &reaction,
                                               const EMU &produced_emu);


std::vector<EMUReaction> CreateAllEMUReactions(const Reaction &reaction,
                                               const EMU &produced_emu);

EMUReaction CreateOneEMUReaction(const Reaction &reaction,
                                 const Substrate &produced_emu_substrate,
                                 const EMU &produced_emu);

std::vector<EMUReaction> SelectUniqueEMUReactions(const std::vector<EMUReaction> &all_new_emu_reactions);

void AddNewEMUInQueue(std::queue<EMU> *queue,
                      const std::set<EMU> &emu_ignore_list,
                      const EMUReactionSide &reaction_side);

#endif //CFLEX_CREATE_EMU_REACTIONS_H

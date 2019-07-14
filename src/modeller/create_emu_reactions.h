#pragma once

#include <vector>
#include <string>
#include <queue>
#include <set>

#include "utilities/reaction.h"
#include "utilities/emu.h"


std::vector<EmuReaction> CreateAllEMUReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<Emu> &observable_emus);


std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const Emu &emu);


std::vector<EmuReaction> CreateNewEMUReactions(const Reaction &reaction,
                                               const Emu &produced_emu);


std::vector<EmuReaction> CreateAllEMUReactions(const Reaction &reaction,
                                               const Emu &produced_emu);

EmuReaction CreateOneEMUReaction(const Reaction &reaction,
                                 const Substrate &produced_emu_substrate,
                                 const Emu &produced_emu);

std::vector<EmuReaction> SelectUniqueEMUReactions(const std::vector<EmuReaction> &all_new_emu_reactions);

void AddNewEMUInQueue(std::queue<Emu> *queue,
                      const std::set<Emu> &emu_ignore_list,
                      const EmuReactionSide &reaction_side);
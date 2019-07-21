#pragma once

#include <vector>
#include <string>
#include <queue>
#include <set>

#include "utilities/reaction.h"
#include "utilities/emu.h"


namespace khnum {
namespace modelling_utills {
// Creates all EMU reactions required for calculation MID of the measured_isotopes
std::vector<EmuReaction> CreateAllEmuReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<Emu> &measured_isotopes);

// Returns every reaction that produce the emu
std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const Emu &emu);

// Creates set of Emu reactions which are produce the emu
std::vector<EmuReaction> CreateNewEmuReactions(const Reaction &reaction,
                                               const Emu &emu);

EmuReaction CreateOneEmuReaction(const Reaction &reaction,
                                 const Substrate &substrate,
                                 const Emu &emu);

std::vector<EmuReaction> SelectUniqueEmuReactions(const std::vector<EmuReaction> &all_new_emu_reactions);

void AddNewEmuInQueue(std::queue<Emu> *queue,
                      const std::set<Emu, decltype(&comparator)> &emu_ignore_list,
                      const EmuReactionSide &reaction_side);
} // namespace modelling_utills
} // namespace khnum
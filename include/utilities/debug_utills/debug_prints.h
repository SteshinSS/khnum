#pragma once

#include <utilities/reaction.h>
#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"

namespace khnum {
    void PrintEmu(const Emu& emu);

    void PrintEmuSubstrate(const EmuSubstrate& emu);

    void PrintEmuReaction(const EmuReaction& reaction);

    void PrintEmuAndMid(const EmuAndMid& emu_and_mid);

    void PrintVectorOfEmu(const std::vector<Emu>& vec);

    void PrintReactions(const std::vector<Reaction>& reactions);

    void PrintReaction(const Reaction& reaction);
} // namespace khnum
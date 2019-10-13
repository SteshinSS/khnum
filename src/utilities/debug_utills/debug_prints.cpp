#include "utilities/debug_utills/debug_prints.h"

#include <iostream>
#include <cmath>

#include "utilities/emu.h"

namespace khnum {
void PrintEmuReaction(const EmuReaction &reaction) {
    std::cout << "id: " << reaction.id << " ";
    bool is_first = true;
    for (const EmuSubstrate& emu : reaction.left) {
        if (!is_first) {
            std::cout << "+ ";
        }
        PrintEmuSubstrate(emu);
        std::cout << " ";
        is_first = false;
    }
    std::cout << "= ";
    PrintEmuSubstrate(reaction.right);
    std::cout << std::endl;
}

void PrintEmuSubstrate(const EmuSubstrate& emu) {
    if (fabs(emu.coefficient - 1.0) > 0.0001) {
        std::cout << emu.coefficient << " ";
    }
    PrintEmu(emu.emu);
}


void PrintEmu(const Emu& emu) {
    std::cout << emu.name << ":";
    for (int atom_state : emu.atom_states) {
        std::cout << atom_state;
    }
}


void PrintEmuAndMid(const EmuAndMid& emu_and_mid) {
    PrintEmu(emu_and_mid.emu);
    std::cout << " [";
    for (double fraction : emu_and_mid.mid) {
        std::cout << fraction << ", ";
    }
    std::cout << "]" << std::endl;
}

void PrintVectorOfEmu(const std::vector<Emu>& vec) {
    for (const Emu& emu : vec) {
        PrintEmu(emu);
        std::cout << std::endl;
    }
}

} // namespace khnum
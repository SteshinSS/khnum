#include "debug_utilites.h"
#include <iostream>
#include <string>

void ShowEMU(const Emu &emu, const std::string &what) {
    std::cerr << emu.name << ":";
    for (bool const &state : emu.atom_states) {
        std::cerr << static_cast<int>(state);
    }
    std::cerr << what;
}


void ShowEMUReaction(const EMUReaction &emu_reaction, const std::string &what) {
    std::cerr << emu_reaction.id << std::endl;
    for (int i = 0; i + 1 < emu_reaction.left.size(); ++i) {
        std::cerr << emu_reaction.left[i].coefficient << " ";
        ShowEMU(emu_reaction.left[i].emu);
        std::cerr << " + ";
    }
    std::cerr << emu_reaction.left.back().coefficient << " ";
    ShowEMU(emu_reaction.left.back().emu);
    std::cerr << " = ";

    std::cerr << emu_reaction.right.coefficient << " ";
    ShowEMU(emu_reaction.right.emu, "\n");

    std::cerr << what;
}
#include "create_emu_list.h"

#include "EMU.h"
#include "input_substrate.h"

#include <vector>
#include <string>

std::vector<EMU> CreateInputEMUList(const std::vector<EMUReaction> &reactions,
                                    const std::vector<InputSubstrate> &input_substrates) {
    std::vector<EMU> input_emu_list;
    for (const EMUReaction &reaction : reactions) {
        // check the left side
        for (const EMUSubstrate &emu_substrate : reaction.left) {
            std::string emu_name = emu_substrate.emu.name;
            for (const InputSubstrate &input_substrate : input_substrates) {
                if (input_substrate.name == emu_name) {
                    input_emu_list.push_back(emu_substrate.emu);
                    break;
                }
            }
        }

        // check the right side
        std::string emu_name = reaction.right.emu.name;
        for (const InputSubstrate &input_substrate : input_substrates) {
            if (input_substrate.name == emu_name) {
                input_emu_list.push_back(reaction.right.emu);
            }
        }
    }
    return input_emu_list;
}
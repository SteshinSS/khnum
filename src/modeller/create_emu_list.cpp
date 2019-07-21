#include "create_emu_list.h"

#include <vector>
#include <string>

#include "utilities/emu.h"
#include "utilities/input_substrate.h"


namespace khnum {
namespace modelling_utills {
std::vector<Emu> CreateInputEmuList(const std::vector<EmuReaction> &reactions,
                                    const std::vector<InputSubstrate> &input_substrates) {
    std::vector<Emu> input_emu_list;
    for (const EmuReaction &reaction : reactions) {
        // check the left side
        for (const EmuSubstrate &emu_substrate : reaction.left) {
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
} // namespace modelling_utills
} // namespace khnum
#include "calculate_input_mid.h"

#include <vector>
#include <algorithm>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/input_substrate.h"


namespace khnum {
namespace modelling_utills {
std::vector<EmuAndMid> CalculateInputMid(const std::vector<InputSubstrate> &input_substrates,
                                         const std::vector<Emu> &input_emus) {
    std::vector<EmuAndMid> input_mids;
    for (const Emu &input_emu : input_emus) {
        // find input substrate with such name
        auto input_substrate_iterator = std::find_if(input_substrates.begin(),
                                                     input_substrates.end(),
                                                     [input_emu](const InputSubstrate &input_substrate) {
                                                         return input_substrate.name == input_emu.name;
                                                     });

        EmuAndMid new_mid = CalculateOneMid(*input_substrate_iterator, input_emu);
        input_mids.push_back(new_mid);
    }
    return input_mids;
}


EmuAndMid CalculateOneMid(const InputSubstrate &input_substrate,
                          const Emu &input_emu) {
    EmuAndMid new_emu_mid;
    new_emu_mid.emu = input_emu;

    // Let input_emu be PYR:1010
    // It means that Emu consists of two traced atoms of PYR: at 0 and 2 positions

    // Find atom positions included in input_emu
    // For our example included_atoms = [0, 2]
    std::vector<int> included_atoms;
    for (size_t i = 0; i < input_emu.atom_states.size(); ++i) {
        if (input_emu.atom_states[i]) {
            included_atoms.push_back(i);
        }
    }

    int emu_size = included_atoms.size();
    // Mid is a vector which i'th value equal fraction of such Emu with mass shift = i
    // For our example there are 3 posibilities of mass shifts: M + 0, M + 1, or M + 2
    Mid new_mid(emu_size + 1, 0.0);
    for (int mass_shift = 0; mass_shift < emu_size + 1; ++mass_shift) {
        double fraction = 0.0;
        for (const auto &mixture : input_substrate.mixtures) {
            // fraction of this input mixture of such mass shift
            double mixture_fraction = 0.0;

            // create bitmask with possible isotope atoms
            // for our example, when mass_shift = 1, we could choose either first or second atoms of PYR:1010
            // so chosen_positions is either [0,1] or [1,0] as well
            std::vector<bool> chosen_positions(emu_size, false);
            for (int position = 0; position < mass_shift; ++position) {
                chosen_positions[emu_size - position - 1] = true;
            }

            // permutate through the all possible masks
            do {
                double current_position_fraction = 1;
                for (int atom_position = 0; atom_position < emu_size; ++atom_position) {
                    if (chosen_positions[atom_position]) {
                        current_position_fraction *= mixture.fractions[included_atoms[atom_position]];
                    } else {
                        current_position_fraction *= 1 - mixture.fractions[included_atoms[atom_position]];
                    }
                }
                mixture_fraction += current_position_fraction;
            } while (std::next_permutation(chosen_positions.begin(), chosen_positions.end()));
            fraction += mixture_fraction * mixture.ratio;
        }
        new_mid[mass_shift] = fraction;
    }
    new_emu_mid.mid = new_mid;

    return new_emu_mid;
}
} // namespace modelling_utills
} // namespace khnum
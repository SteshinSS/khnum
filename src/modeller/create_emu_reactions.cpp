#include "modeller/create_emu_reactions.h"

#include <vector>
#include <string>
#include <queue>
#include <exception>
#include <set>
#include <algorithm>
#include <iostream>

#include "utilities/emu.h"
#include "utilities/reaction.h"
#include "utilities/debug_utills/debug_prints.h"


namespace khnum {
namespace modelling_utills {
std::vector<EmuReaction> CreateAllEmuReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<Emu> &measured_isotopes) {
    std::vector<EmuReaction> all_emu_reactions;

    // contains all EMUs for finding their synthesis reactions
    std::queue<Emu> emu_to_check;

    // contains checked EMUs
    std::set<Emu, decltype(&comparator)> already_checked_emu(&comparator);

    for (const Emu &emu : measured_isotopes) { // initializing queue
        emu_to_check.push(emu);
    }

    while (!emu_to_check.empty()) {
        Emu next_emu = emu_to_check.front();
        emu_to_check.pop();

        if (already_checked_emu.find(next_emu) == already_checked_emu.end()) {
            // find all reactions which synthesize the next_emu
            std::vector<Reaction> synthesis_reactions = GetSynthesisReactions(reactions, next_emu);

            for (const Reaction &reaction : synthesis_reactions) {
                std::vector<EmuReaction> new_emu_reactions = CreateNewEmuReactions(reaction, next_emu);
                for (const EmuReaction &emu_reaction : new_emu_reactions) {
                    all_emu_reactions.push_back(emu_reaction);
                    AddNewEmusInQueue(&emu_to_check, already_checked_emu, emu_reaction.left);
                }
            }
            already_checked_emu.insert(next_emu);
        }
    }
    return all_emu_reactions;
}


std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const Emu &emu) {
    std::vector<Reaction> synthesis_reactions;
    for (const Reaction &reaction : reactions) {
        if (reaction.type != ReactionType::MetaboliteBalance) {
            bool is_reaction_produce_emu = false;
            for (const Substrate &substrate : reaction.chemical_equation.right) {
                if (substrate.name == emu.name) {
                    is_reaction_produce_emu = true;
                    break;
                }
            }

            if (is_reaction_produce_emu) {
                synthesis_reactions.push_back(reaction);
            }
        }
    }
    return synthesis_reactions;
}


std::vector<EmuReaction> CreateNewEmuReactions(const Reaction &reaction,
                                               const Emu &emu) {
    std::vector<EmuReaction> new_emu_reactions;

    // find all occurrences of the emu in the reaction
    for (const Substrate &substrate : reaction.chemical_equation.right) {
        if (substrate.name == emu.name) {
            new_emu_reactions.push_back(CreateOneEmuReaction(reaction, substrate, emu));
        }
    }

    std::vector<EmuReaction> unique_emu_reactions = SelectUniqueEmuReactions(new_emu_reactions);
    return unique_emu_reactions;
}


EmuReaction CreateOneEmuReaction(const Reaction &reaction,
                                 const Substrate &substrate,
                                 const Emu &emu) {
    EmuReaction result_reaction;
    result_reaction.id = reaction.id;

    // form right side
    EmuSubstrate product;
    product.emu = emu;
    product.coefficient = substrate.coefficient;
    result_reaction.right = product;

    EmuReactionSide left_side;

    // contains all substrates that already included in the equation's left side
    // and their position
    std::vector<std::pair<Substrate, int>> included_substrates;

    // form left side
    // find sources of all atoms included in the emu
    for (size_t atom_position = 0; atom_position < emu.atom_states.size(); ++atom_position) {
        if (emu.atom_states[atom_position]) {
            char atom = substrate.formula[atom_position]; // letter represents atom in the formula

            // finding precursor substrate that produces this atom
            bool is_atom_found = false;
            for (const Substrate &precursor : reaction.chemical_equation.left) {
                size_t substrate_atom_position = precursor.formula.find(atom);
                if (substrate_atom_position != std::string::npos) {
                    is_atom_found = true;

                    // checking if this precursor substrate is already in our elementary reaction
                    auto precursor_iterator = std::find_if(included_substrates.begin(),
                                                           included_substrates.end(),
                                                           [&precursor](const std::pair<Substrate, int> pair) {
                                                               return pair.first == precursor;
                                                           });

                    if (precursor_iterator != included_substrates.end()) {
                        int precursor_position = precursor_iterator->second;
                        left_side[precursor_position].emu.atom_states[substrate_atom_position] = true;
                    } else {
                        EmuSubstrate new_precursor;
                        new_precursor.emu.name = precursor.name;
                        new_precursor.emu.atom_states = AtomStates(precursor.formula.size(), false);
                        new_precursor.emu.atom_states[substrate_atom_position] = true;
                        new_precursor.coefficient = precursor.coefficient;
                        left_side.push_back(new_precursor);
                        included_substrates.push_back({precursor, left_side.size() - 1});
                    }
                    break;
                }
            }

            if (!is_atom_found) {
                throw std::runtime_error("There is reaction in which there are atoms "
                                         "in the right side such that are not in the left one " +
                    std::to_string(reaction.id));

            }
        }
    }
    result_reaction.left = left_side;
    return result_reaction;
}


std::vector<EmuReaction> SelectUniqueEmuReactions(const std::vector<EmuReaction> &emu_reactions) {
    std::vector<EmuReaction> unique_emu_reactions;
    for (const EmuReaction &reaction : emu_reactions) {
        auto reaction_position = std::find(unique_emu_reactions.begin(),
                                           unique_emu_reactions.end(), reaction);
        if (reaction_position == unique_emu_reactions.end()) {
            unique_emu_reactions.push_back(reaction);
        } else {
            (reaction_position->right).coefficient += reaction.right.coefficient;
        }
    }

    return unique_emu_reactions;
}


void AddNewEmusInQueue(std::queue<Emu> *queue,
                       const std::set<Emu, decltype(&comparator)> &already_checked_emu,
                       const EmuReactionSide &reaction_side) {
    for (EmuSubstrate const &substrate : reaction_side) {
        if (already_checked_emu.find(substrate.emu) == already_checked_emu.end()) {
            queue->push(substrate.emu);
        }
    }

}
} // namespace modelling_utills
} // namespace khnum
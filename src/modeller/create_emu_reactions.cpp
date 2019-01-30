#include "create_emu_reactions.h"
#include "../utilities/EMU.h"
#include "../utilities/reaction_struct.h"
#include "../utilities/debug_utilites.h"

#include <vector>
#include <string>
#include <queue>
#include <exception>
#include <set>
#include <algorithm>
#include <iostream>
std::vector<EMUReaction> CreateAllEMUReactions(const std::vector<Reaction> &reactions,
                                               const std::vector<EMU> &observable_emus) {
    std::vector<EMUReaction> all_emu_reactions;
    std::queue<EMU> emu_to_check; // contains all EMUs for finding their synthesis reactions
    std::set<EMU> already_checked_emu; // contains checked EMUs
    for (const EMU &emu : observable_emus) { // initializing queue
        emu_to_check.push(emu);
    }

    while (!emu_to_check.empty()) {
        EMU next_emu = emu_to_check.front();
        emu_to_check.pop();

        if (already_checked_emu.find(next_emu) == already_checked_emu.end()) {

            // find all reactions which synthesize the next_emu
            std::vector<Reaction> synthesis_reactions = GetSynthesisReactions(reactions, next_emu);

            for (const Reaction &reaction : synthesis_reactions) {
                std::vector<EMUReaction> new_emu_reactions = CreateNewEMUReactions(reaction, next_emu);
                for (const EMUReaction &emu_reaction : new_emu_reactions) {
                    all_emu_reactions.push_back(emu_reaction);
                    AddNewEMUInQueue(&emu_to_check, already_checked_emu, emu_reaction.left);
                }
            }
            already_checked_emu.insert(next_emu);
        }
    }
    return all_emu_reactions;
}

std::vector<Reaction> GetSynthesisReactions(const std::vector<Reaction> &reactions,
                                            const EMU &emu) {
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

// Creates set of EMU reactions which are produce the produced_emu
std::vector<EMUReaction> CreateNewEMUReactions(const Reaction &reaction,
                                               const EMU &produced_emu) {
    std::vector<EMUReaction> new_emu_reactions = CreateAllEMUReactions(reaction, produced_emu);
    std::vector<EMUReaction> unique_emu_reactions = SelectUniqueEMUReactions(new_emu_reactions);
    return unique_emu_reactions;
}

// Creates set of all EMU reactions. Repetitions are possible
std::vector<EMUReaction> CreateAllEMUReactions(const Reaction &reaction,
                                               const EMU &produced_emu) {
    std::vector<EMUReaction> all_new_emu_reactions;
    // find all occurrences of produced_emu in the reaction
    for (const Substrate &substrate : reaction.chemical_equation.right) {
        if (substrate.name == produced_emu.name) {
            all_new_emu_reactions.push_back(CreateOneEMUReaction(reaction, substrate, produced_emu));
        }
    }
    return all_new_emu_reactions;
}

EMUReaction CreateOneEMUReaction(const Reaction &reaction,
                                 const Substrate &produced_emu_substrate,
                                 const EMU &produced_emu) {
    EMUReaction result_reaction;
    result_reaction.id = reaction.id;

    // form right side
    EMUSubstrate product;
    product.emu = produced_emu;
    product.coefficient = produced_emu_substrate.coefficient;
    result_reaction.right = product;

    // contain all substrates that already included in the equation's left side
    // and their position
    std::vector<std::pair<Substrate, int>> included_substrates;

    // form left side
    // find all atoms included in produced_emu
    EMUReactionSide left;
    for (int atom_position = 0; atom_position < produced_emu.atom_states.size(); ++atom_position) {
        if (produced_emu.atom_states[atom_position]) {

            char atom_name = produced_emu_substrate.formula[atom_position]; // letter representing atom in formula

            // finding precursor substrate that produces this atom
            bool is_atom_found = false;
            for (const Substrate &precursor : reaction.chemical_equation.left) {
                int substrate_atom_position = precursor.formula.find(atom_name);
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
                        left[precursor_position].emu.atom_states[substrate_atom_position] = true;
                    } else {
                        EMUSubstrate new_precursor;
                        new_precursor.emu.name = precursor.name;
                        new_precursor.emu.atom_states = AtomStates(precursor.formula.size(), false);
                        new_precursor.emu.atom_states[substrate_atom_position] = true;
                        new_precursor.coefficient = precursor.coefficient;
                        left.push_back(new_precursor);
                        included_substrates.push_back({precursor, left.size() - 1});
                    }
                    break;
                }
            }

            if (!is_atom_found) {
                throw std::runtime_error("There is reaction in which there are atoms "
                                         "in the right side such that are not in the left one");
            }
        }
    }

    result_reaction.left = left;
    return result_reaction;
}

std::vector<EMUReaction> SelectUniqueEMUReactions(const std::vector<EMUReaction> &all_new_emu_reactions) {
    std::vector<EMUReaction> unique_new_emu_reactions;
    for (const EMUReaction &reaction : all_new_emu_reactions) {
        auto reaction_position = std::find(unique_new_emu_reactions.begin(),
                                           unique_new_emu_reactions.end(), reaction);
        if (reaction_position == unique_new_emu_reactions.end()) {
            unique_new_emu_reactions.push_back(reaction);
        } else {
            (reaction_position->right).coefficient += reaction.right.coefficient;
        }
    }

    return unique_new_emu_reactions;
}

void AddNewEMUInQueue(std::queue<EMU> *queue,
                      const std::set<EMU> &emu_ignore_list,
                      const EMUReactionSide &reaction_side) {
    for (EMUSubstrate const &substrate : reaction_side) {
        if (emu_ignore_list.find(substrate.emu) == emu_ignore_list.end()) {
            queue->push(substrate.emu);
        }
    }

}
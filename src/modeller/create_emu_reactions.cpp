#include "modeller/create_emu_reactions.h"

#include <vector>
#include <string>
#include <queue>
#include <exception>
#include <set>
#include <algorithm>
#include <iostream>
#include <map>

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
    Rate rate = 1.0;

    auto substrate_iterator = std::find(reaction.chemical_equation.right.begin(),
                                        reaction.chemical_equation.right.end(), substrate);

    const int substrate_pos = substrate_iterator - reaction.chemical_equation.right.begin();

    // form right side
    EmuSubstrate product;
    product.emu = emu;
    product.coefficient = substrate.substrate_coefficient_;
    result_reaction.right = product;
    rate *= product.coefficient;

    EmuReactionSide left_side;

    // contains all substrates that already included in the equation's left side
    // and their position
    std::map<int, int> reactionpos2emupos;
    // form left side
    // find sources of all atoms included in the emu
    for (size_t atom_position = 0; atom_position < emu.atom_states.size(); ++atom_position) {
        if (emu.atom_states[atom_position]) {
            std::vector<AtomTransition> transitions;
            for (const AtomTransition &transition: reaction.chemical_equation.atom_transitions) {
                if (transition.product_pos == substrate_pos && transition.product_atom == atom_position) {
                    transitions.push_back(transition);
                }
            }

            for (const AtomTransition &transition : transitions) {
                if (reactionpos2emupos.find(transition.substrate_pos) != reactionpos2emupos.end()) {
                    left_side[reactionpos2emupos[transition.substrate_pos]].emu.atom_states[transition.substrate_atom] = true;
                } else {
                    const Substrate &precursor = reaction.chemical_equation.left[transition.substrate_pos];
                    EmuSubstrate new_precursor;
                    new_precursor.emu.name = precursor.name;
                    new_precursor.emu.atom_states = AtomStates(precursor.size, false);
                    new_precursor.emu.atom_states[transition.substrate_atom] = true;
                    new_precursor.coefficient = precursor.substrate_coefficient_;
                    left_side.push_back(new_precursor);
                    reactionpos2emupos.insert({transition.substrate_pos, left_side.size() - 1});
                }
            }
        }
    }
    for (const EmuSubstrate& precursor : result_reaction.left) {
        rate *= precursor.coefficient;
    }

    result_reaction.left = left_side;
    result_reaction.rate = rate;

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
            reaction_position->rate += reaction.rate;
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

#include "sort_reactions.h"


#include "../utilities/reaction_struct.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

#include <iostream>


std::vector<Reaction> SortReactions(std::vector<Reaction> reactions) {
    std::sort(reactions.begin(), reactions.end(), [](const Reaction &lhs, const Reaction &rhs) {
        return GetPriority(lhs) < GetPriority(rhs);
    });

    return reactions;
}

int GetPriority(const Reaction &reaction) {
    int priority = -1;
    if (std::isnan(reaction.basis) && !reaction.is_basis_x) {
        if (reaction.type == ReactionType::Forward) {
            priority = 1;
        } else if (reaction.type == ReactionType::Irreversible) {
            priority = 2;
        } else if (reaction.type == ReactionType::Backward) {
            priority = 6;
        } else if (reaction.type == ReactionType::MetaboliteBalance) {
            priority = 3;
        } else if (reaction.type == ReactionType::IsotopomerBalance) {
            priority = 0;
        }
    } else if (reaction.is_basis_x) {
        priority = 4;
    } else {
        priority = 5;
    }

    return priority;
}



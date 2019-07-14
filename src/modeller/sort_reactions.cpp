#include "sort_reactions.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

#include "utilities/reaction.h"


std::vector<Reaction> SortReactionsByType(std::vector<Reaction> reactions) {
    std::sort(reactions.begin(), reactions.end(), [](const Reaction &lhs, const Reaction &rhs) {
        int lhs_priority = GetPriority(lhs);
        int rhs_priority = GetPriority(rhs);
        if (lhs_priority == rhs_priority) {
            return lhs.id < rhs.id;
        } else {
            return lhs_priority < rhs_priority;
        }
    });

    return reactions;
}

int GetPriority(const Reaction &reaction) {
    int priority = -1;
    if (std::isnan(reaction.basis) && !reaction.is_set_free) {
        if (reaction.type == ReactionType::Forward) {
            priority = 1;
        } else if (reaction.type == ReactionType::Irreversible) {
            priority = 2;
        } else if (reaction.type == ReactionType::Backward) {
            priority = 4;
        } else if (reaction.type == ReactionType::MetaboliteBalance) {
            priority = 3;
        } else if (reaction.type == ReactionType::IsotopomerBalance) {
            priority = 0;
        }
    } else if (reaction.is_set_free) {
        priority = 5;
    } else {
        priority = 6;
    }

    return priority;
}


std::vector<Reaction> SortReactionByID(std::vector<Reaction> reactions) {
    std::sort(reactions.begin(), reactions.end(), [](const Reaction &lhs, const Reaction &rhs) {
        return lhs.id < rhs.id;
    });

    return reactions;
}


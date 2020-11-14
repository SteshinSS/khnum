#include "modeller/sort_reactions.h"

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>

#include "utilities/reaction.h"


namespace khnum {
namespace modelling_utills {
std::vector<Reaction> SortReactionsByType(const std::vector<Reaction>& reactions) {
    std::vector<Reaction> sorted_reactions = reactions;
    std::sort(sorted_reactions.begin(), sorted_reactions.end(), [](const Reaction &lhs, const Reaction &rhs) {
        int lhs_priority = GetPriority(lhs);
        int rhs_priority = GetPriority(rhs);
        if (lhs_priority == rhs_priority) {
            return lhs.id > rhs.id;
        } else {
            return lhs_priority < rhs_priority;
        }
    });

    return sorted_reactions;
}


int GetPriority(const Reaction &reaction) {
    int priority = -1;
    if (std::isnan(reaction.basis) && !reaction.is_set_free) {
        if (reaction.type == ReactionType::Forward) {
            priority = 1;
        } else if (reaction.type == ReactionType::Irreversible) {
            priority = 3;
        } else if (reaction.type == ReactionType::Backward) {
            priority = 4;
        } else if (reaction.type == ReactionType::MetaboliteBalance) {
            priority = 2;
        } else if (reaction.type == ReactionType::IsotopomerBalance) {
            priority = 0;
        }
    } else if (std::isnan(reaction.basis)) {
        priority = 6;
    } else {
        priority = 5;
    }

    return priority;
}


std::vector<Reaction> SortReactionByID(std::vector<Reaction> reactions) {
    std::sort(reactions.begin(), reactions.end(), [](const Reaction &lhs, const Reaction &rhs) {
        return lhs.id < rhs.id;
    });

    return reactions;
}
} // namespace modelling_utills
} // namespace khnum
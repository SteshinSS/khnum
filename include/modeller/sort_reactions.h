#pragma once

#include "utilities/reaction.h"

#include <vector>
#include <algorithm>
#include <limits>


namespace khnum {
namespace modelling_utills {
std::vector<Reaction> SortReactionsByType(std::vector<Reaction> reactions);

std::vector<Reaction> SortReactionByID(std::vector<Reaction> reactions);

int GetPriority(const Reaction &reaction);
} // namespace modelling_utills
} // namespace khnum
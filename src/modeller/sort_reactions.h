#pragma once

#include "utilities/reaction.h"

#include <vector>
#include <algorithm>
#include <limits>

std::vector<Reaction> SortReactionsByType(std::vector<Reaction> reactions);

std::vector<Reaction> SortReactionByID(std::vector<Reaction> reactions);

int GetPriority(const Reaction &reaction);
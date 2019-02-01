#ifndef CFLEX_SORT_REACTIONS_H
#define CFLEX_SORT_REACTIONS_H


#include "reaction_struct.h"

#include <vector>
#include <algorithm>
#include <limits>

std::vector<Reaction> SortReactionsByType(std::vector<Reaction> reactions);

std::vector<Reaction> SortReactionByID(std::vector<Reaction> reactions);

int GetPriority(const Reaction &reaction);

#endif //CFLEX_SORT_REACTIONS_H

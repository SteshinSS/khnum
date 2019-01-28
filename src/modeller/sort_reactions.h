#ifndef CFLEX_SORT_REACTIONS_H
#define CFLEX_SORT_REACTIONS_H


#include "../utilities/reaction_struct.h"

#include <vector>
#include <algorithm>
#include <limits>

std::vector<Reaction> SortReactions(std::vector<Reaction> reactions);


int GetPriority(const Reaction &reaction);

#endif //CFLEX_SORT_REACTIONS_H

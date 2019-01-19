#ifndef CFLEX_CREATE_METABOLITE_LIST_H
#define CFLEX_CREATE_METABOLITE_LIST_H

#include "../utilities/reaction_struct.h"

#include <vector>
#include <string>
#include <algorithm>

std::vector<std::string> CreateFullMetaboliteList(const std::vector<Reaction> &reactions);

std::vector<std::string> CreateIncludedMetaboliteList(const std::vector<std::string> &metabolite_list,
                                                      const std::vector<std::string> &excluded_metabolites);


#endif //CFLEX_CREATE_METABOLITE_LIST_H

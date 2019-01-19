#include "create_metabolite_list.h"

#include <vector>
#include <string>
#include <algorithm>

std::vector<std::string> CreateFullMetaboliteList(const std::vector<Reaction> &reactions) {
    std::vector<std::string> metabolite_list;
    for (auto const &reaction : reactions) {
        if (reaction.type != ReactionType::IsotopomerBalance) {
            for (auto const &substance : reaction.chemical_equation.left) {
                metabolite_list.push_back(substance.name);
            }
            for (auto const &substance : reaction.chemical_equation.right) {
                metabolite_list.push_back(substance.name);
            }
        }

    }
    std::sort(metabolite_list.begin(), metabolite_list.end());
    metabolite_list.erase(std::unique(metabolite_list.begin(), metabolite_list.end()), metabolite_list.end());
    return metabolite_list;
}

std::vector<std::string> CreateIncludedMetaboliteList(const std::vector<std::string> &metabolite_list,
                                                      const std::vector<std::string> &excluded_metabolites) {
    std::vector<std::string> included_metabolite_list;
    for (auto const &metabolite : metabolite_list) {
        if (std::find(excluded_metabolites.begin(), excluded_metabolites.end(), metabolite) ==
            excluded_metabolites.end()) {
            included_metabolite_list.emplace_back(metabolite);
        }
    }
    return included_metabolite_list;
}
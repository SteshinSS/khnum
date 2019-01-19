#include "create_stoichiometry_matrix.h"
#include "../math/math_utilites.h"
#include "../utilities/reaction_struct.h"

#include <string>
#include <vector>
#include <exception>

Matrix CreateStoichiometryMatrix(const std::vector<Reaction> &reactions,
                                 const std::vector<std::string> &metabolite_list) {

    const int metabolite_number = metabolite_list.size();

    // counting number of reaction included in stoichiometry matrix
    int reaction_number = 0;
    for (int reaction = 0; reaction < reactions.size(); ++reaction) {
        if (reactions[reaction].type != ReactionType::IsotopomerBalance) {
            ++reaction_number;
        }
    }
    Matrix stoichiometry_matrix(metabolite_number, reaction_number);

    for (int reaction = 0; reaction < reaction_number; ++reaction) {
        if (reactions[reaction].type != ReactionType::IsotopomerBalance) {
            for (int metabolite = 0; metabolite < metabolite_number; ++metabolite) {
                std::string current_metabolite_name = metabolite_list.at(metabolite);
                stoichiometry_matrix(metabolite, reaction) = GetTotalCoefficient(
                    reactions.at(reaction).chemical_equation,
                    current_metabolite_name);
            }
        }
    }

    return stoichiometry_matrix;
}

double GetTotalCoefficient(const ChemicalEquation &chemical_equation, const std::string &metabolite) {
    double result{0.0};
    bool is_found_at_left{false};
    for (auto const &substrate : chemical_equation.left) {
        if (substrate.name == metabolite) {
            result -= substrate.coefficient;
            is_found_at_left = true;
        }
    }
    for (auto const &substrate : chemical_equation.right) {
        if (substrate.name == metabolite) {
            if (!is_found_at_left) {
                result += substrate.coefficient;
            } else {
                throw std::runtime_error("There is reaction with the same substrate in both sides!");
            }
        }
    }
    return result;
}
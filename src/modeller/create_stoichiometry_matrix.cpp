#include "create_stoichiometry_matrix.h"
#include "math_utilites.h"
#include "reaction_struct.h"

#include <string>
#include <vector>
#include <iostream>
#include <exception>

Matrix CreateStoichiometryMatrix(const std::vector<Reaction> &reactions,
                                 const std::vector<std::string> &metabolite_list) {

    const int metabolite_number = metabolite_list.size();

    // counting number of reaction included in stoichiometry matrix
    int reaction_number = 0;
    for (const Reaction &reaction : reactions) {
        if (reaction.type != ReactionType::IsotopomerBalance) {
            ++reaction_number;
        }
    }
    Matrix stoichiometry_matrix = Matrix::Zero(metabolite_number, reaction_number);

    int stoichiometry_reaction_number = 0;
    for (const Reaction &reaction : reactions) {
        if (reaction.type != ReactionType::IsotopomerBalance) {
            for (int metabolite = 0; metabolite < metabolite_number; ++metabolite) {
                std::string current_metabolite_name = metabolite_list.at(metabolite);
                stoichiometry_matrix(metabolite, stoichiometry_reaction_number) = GetTotalCoefficient(
                    reaction.chemical_equation,
                    current_metabolite_name);
            }
            ++stoichiometry_reaction_number;
        }
    }

    return stoichiometry_matrix;
}

double GetTotalCoefficient(const ChemicalEquation &chemical_equation, const std::string &metabolite) {
    double result{0.0};
    bool is_found_at_left{false};
    for (const Substrate &substrate : chemical_equation.left) {
        if (substrate.name == metabolite) {
            result -= substrate.coefficient;
            is_found_at_left = true;
        }
    }
    for (const Substrate &substrate : chemical_equation.right) {
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
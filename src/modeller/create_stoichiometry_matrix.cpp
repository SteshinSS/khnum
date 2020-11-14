#include "modeller/create_stoichiometry_matrix.h"

#include <string>
#include <vector>
#include <iostream>
#include <exception>

#include "utilities/matrix.h"
#include "utilities/reaction.h"
#include "utilities/debug_utills/debug_prints.h"


namespace khnum {
namespace modelling_utills {
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
        bool is_present = false;
        if (reaction.type != ReactionType::IsotopomerBalance) {
            for (int metabolite = 0; metabolite < metabolite_number; ++metabolite) {
                std::string current_metabolite_name = metabolite_list.at(metabolite);

                double coefficient = GetTotalCoefficient(reaction.chemical_equation, current_metabolite_name);
                if (coefficient < -0.001 || coefficient > 0.001) {
                    is_present = true;
                }
                stoichiometry_matrix(metabolite, stoichiometry_reaction_number) = coefficient;
            }
            ++stoichiometry_reaction_number;
        }
        if (!is_present) {
            std::cout << "Reaction with no coefficients: " << std::endl;
            PrintReaction(reaction);
        }
    }

    // Check for metabolites which doesn't have input and output together
    for (int i = 0; i < stoichiometry_matrix.rows(); ++i) {
        bool is_pos = false;
        bool is_neg = false;
        for (int j = 0; j < stoichiometry_matrix.cols(); ++j) {
            if (stoichiometry_matrix(i, j) > 0.0) {
                is_pos = true;
            }
            if (stoichiometry_matrix(i, j) < -0.0) {
                is_neg = true;
            }
        }
        if (!(is_neg && is_pos)) {
            std::cout << "bad: " << metabolite_list[i] << std::endl;
        }
    }

    return stoichiometry_matrix;
}


double GetTotalCoefficient(const ChemicalEquation &chemical_equation, const std::string &metabolite) {
    double result{0.0};
    for (const Substrate &substrate : chemical_equation.left) {
        if (substrate.name == metabolite) {
            result -= substrate.substrate_coefficient_;
        }
    }
    for (const Substrate &substrate : chemical_equation.right) {
        if (substrate.name == metabolite) {
            result += substrate.substrate_coefficient_;
        }
    }
    return result;
}
} // namespace modelling_utills
} // namespace khnum
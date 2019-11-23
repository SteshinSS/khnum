#include "modeller/create_stoichiometry_matrix.h"

#include <string>
#include <vector>
#include <iostream>
#include <exception>

#include "utilities/matrix.h"
#include "utilities/reaction.h"


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
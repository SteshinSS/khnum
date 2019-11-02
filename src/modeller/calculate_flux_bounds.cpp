#include "modeller/calculate_flux_bounds.h"

#include <vector>
#include <glpk/glpk.h>

#include "utilities/reaction.h"
#include "utilities/matrix.h"


namespace khnum {
namespace modelling_utills {


// We find upper bound for ith flux in such way:

// maximize w_i * x
// subject to S * v = 0
// where x > 0
// and w_i = (0, 0, ... 0, 1, 0, ... 0) with 1 in ith position

// First reactions in the vector are metabolic balance reactions, so we don't calculate bounds for them
void CalculateFluxBounds(std::vector<Reaction>& reactions, const Matrix& stoichiometry_matrix) {
    glp_term_out(GLP_OFF);
    glp_prob* linear_problem = glp_create_prob();

    PrepareLinearProblem(reactions, stoichiometry_matrix, linear_problem);

    const int total_reactions_with_bounds = stoichiometry_matrix.cols();
    const int metabolite_balance_reactions_total = reactions.size() - total_reactions_with_bounds;
    for (int reaction = 0; reaction < total_reactions_with_bounds; ++reaction) {
        // Set w_i vector
        for (int i = 0; i < total_reactions_with_bounds; ++i) {
            if (i == reaction) {
                glp_set_obj_coef(linear_problem, i + 1, 1.0);
            } else {
                glp_set_obj_coef(linear_problem, i + 1, 0.0);
            }
        }
        glp_set_obj_dir(linear_problem, GLP_MAX);
        glp_simplex(linear_problem, NULL);
        reactions[reaction + metabolite_balance_reactions_total].computed_upper_bound = glp_get_obj_val(linear_problem);

        glp_set_obj_dir(linear_problem, GLP_MIN);
        glp_simplex(linear_problem, NULL);
        reactions[reaction + metabolite_balance_reactions_total].computed_lower_bound = glp_get_obj_val(linear_problem);
    }

    glp_delete_prob(linear_problem);
}

void PrepareLinearProblem(std::vector<Reaction>& reactions, const Matrix& stoichiometry_matrix, glp_prob* linear_problem) {
    // Set LB < xi < UB
    glp_add_cols(linear_problem, stoichiometry_matrix.cols());
    const int total_reactions_with_bounds = stoichiometry_matrix.cols();
    const int metabolite_balance_reactions_total = reactions.size() - total_reactions_with_bounds;
    for (int reaction_num = 0; reaction_num < total_reactions_with_bounds; ++reaction_num) {
        const Reaction& reaction = reactions.at(reaction_num + metabolite_balance_reactions_total);
        if (std::isnan(reaction.basis)) {
            double upper_bound = reaction.setted_upper_bound ?
                                 *reaction.setted_upper_bound : 200;
            double lower_bound = reaction.setted_lower_bound ?
                                 *reaction.setted_lower_bound : 0;
            glp_set_col_bnds(linear_problem, reaction_num + 1, GLP_DB, lower_bound, upper_bound);
        } else {
            if (std::isnan(reaction.deviation)) {
                glp_set_col_bnds(linear_problem, reaction_num + 1, GLP_FX, reaction.basis, reaction.basis);
            } else {
                double upper_bound  = reaction.basis + reaction.deviation;
                double lower_bound = reaction.basis - reaction.deviation;
                glp_set_col_bnds(linear_problem, reaction_num + 1, GLP_DB, lower_bound, upper_bound);
            }
        }
    }

    // Set constraint that S*v = 0
    glp_add_rows(linear_problem, stoichiometry_matrix.rows());
    for (int row = 0; row < stoichiometry_matrix.rows(); ++row) {
        glp_set_row_bnds(linear_problem, row + 1, GLP_FX, 0.0, 0.0);
    }

    // Fill S matrix
    std::vector<int> row_index = {0};
    std::vector<int> col_index = {0};
    std::vector<double> coefficients = {0.0};
    for (int row = 0; row < stoichiometry_matrix.rows(); ++row) {
        for (int col = 0; col < stoichiometry_matrix.cols(); ++col) {
            row_index.push_back(row + 1);
            col_index.push_back(col + 1);
            coefficients.push_back(stoichiometry_matrix(row, col));
        }
    }
    glp_load_matrix(linear_problem, (coefficients.size() - 1), row_index.data(),
                    col_index.data(),
                    coefficients.data());
}
}
}
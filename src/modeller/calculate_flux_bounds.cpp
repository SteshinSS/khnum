#include "modeller/calculate_flux_bounds.h"

#include <vector>
#include <glpk/glpk.h>

#include "utilities/reaction.h"
#include "utilities/matrix.h"


namespace khnum {
namespace modelling_utills {
void CalculateFluxBounds(std::vector<Reaction>& reactions, const Matrix& stoichiometry_matrix) {
    glp_prob* linear_problem = glp_create_prob();
    glp_term_out(GLP_OFF);
    glp_add_rows(linear_problem, stoichiometry_matrix.rows());
    for (int row = 0; row < stoichiometry_matrix.rows(); ++row) {
        glp_set_row_bnds(linear_problem, row + 1, GLP_FX, 0.0, 0.0);
    }

    int metabolite_balance_reactions_total = reactions.size() - stoichiometry_matrix.cols();
    glp_add_cols(linear_problem, stoichiometry_matrix.cols());
    for (int column = 0; column < stoichiometry_matrix.cols(); ++column) {
        const Reaction& reaction = reactions.at(column + metabolite_balance_reactions_total);
        if (std::isnan(reaction.basis)) {
            double upper_bound = reaction.setted_upper_bound ?
                                 *reaction.setted_upper_bound : 10;
            double lower_bound = reaction.setted_lower_bound ?
                                 *reaction.setted_lower_bound : 0;
            glp_set_col_bnds(linear_problem, column + 1, GLP_DB, lower_bound, upper_bound);
        } else {
            if (std::isnan(reaction.deviation)) {
                glp_set_col_bnds(linear_problem, column + 1, GLP_FX, reaction.basis, reaction.basis);
            } else {
                double upper_bound  = reaction.basis + reaction.deviation;
                double lower_bound = reaction.basis - reaction.deviation;
                glp_set_col_bnds(linear_problem, column + 1, GLP_DB, lower_bound, upper_bound);
            }
        }

    }

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
    for (int reaction = 0; reaction < stoichiometry_matrix.cols(); ++reaction) {
        for (int column = 0; column < stoichiometry_matrix.cols(); ++column) {
            if (column == reaction) {
                glp_set_obj_coef(linear_problem, column + 1, 1.0);
            } else {
                glp_set_obj_coef(linear_problem, column + 1, 0.0);
            }
        }
        glp_set_obj_dir(linear_problem, GLP_MAX);
        glp_simplex(linear_problem, NULL);
        reactions[reaction + metabolite_balance_reactions_total].computed_upper_bound = glp_get_obj_val(linear_problem);
        glp_set_obj_dir(linear_problem, GLP_MIN);
        glp_simplex(linear_problem, NULL);
        reactions[reaction + metabolite_balance_reactions_total].computed_lower_bound = glp_get_obj_val(linear_problem);
    }
}
}
}
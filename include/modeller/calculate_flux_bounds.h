#pragma once

#include <vector>
#include <glpk/glpk.h>

#include "utilities/reaction.h"
#include "utilities/matrix.h"

namespace khnum {
namespace modelling_utills {

void CalculateFluxBounds(std::vector<Reaction>& reactions, const Matrix& stoichiometry_matrix);
void PrepareLinearProblem(std::vector<Reaction>& reactions, const Matrix& stoichiometry_matrix, glp_prob* linear_problem);

}
}
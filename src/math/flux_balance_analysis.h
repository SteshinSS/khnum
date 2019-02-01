#ifndef CFLEX_FIND_INITIAL_FLUXES_H
#define CFLEX_FIND_INITIAL_FLUXES_H

#include "reaction_struct.h"
#include "math_utilites.h"
#include "linear_problem.h"

#include <vector>
#include <string>
#include <map>

// Contains FBA-calculated flux ranges
// probably should rename, so it won't confuse with input bounds

void CreateLinearProblem(LinearProblem &linear_problem,
                         const Matrix &stoichiometry_matrix,
                         const std::vector<Reaction> &reactions,
                         const std::vector<std::string> &included_metabolites);

void PrepareMatrixForGLPK(const Matrix &matrix,
                          LinearProblem &linear_problem);

#endif //CFLEX_FIND_INITIAL_FLUXES_H

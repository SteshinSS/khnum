//
// Created by Олег on 2019-01-03.
//

#ifndef CFLEX_FIND_INITIAL_FLUXES_H
#define CFLEX_FIND_INITIAL_FLUXES_H

#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"
#include "../utilities/linear_problem.h"

#include <vector>
#include <string>
#include <map>

// Contains FBA-calculated flux ranges
// probably should rename, so it won't confuse with input bounds
struct FluxVariability {
  Flux lower_bound;
  Flux upper_bound;
};

std::map<std::string, Flux> EstablishInitialFluxes(const Matrix &stoichiometry_matrix,
                                         const std::vector<Reaction> &reactions,
                                         const std::vector<std::string> &included_metabolites);

std::map<std::string, FluxVariability> EstablishAllFluxRanges(const Matrix &stoichiometry_matrix,
                                                    const std::vector<Reaction> &reactions,
                                                    const std::vector<std::string> &included_metabolites);


FluxVariability EstablishFluxRange(int reaction_index,
                                   const Matrix &stoichiometry_matrix,
                                   const std::vector<Reaction> &reactions,
                                   const std::vector<std::string> &included_metabolites);

Flux EstablishExtremeFlux(int reaction_index,
                          const Matrix &stoichiometry_matrix,
                          const std::vector<Reaction> &reactions,
                          const std::vector<std::string> &included_metabolites,
                          bool maximize);

void CreateLinearProblem(LinearProblem &linear_problem,
                         const Matrix &stoichiometry_matrix,
                         const std::vector<Reaction> &reactions,
                         const std::vector<std::string> &included_metabolites);

void PrepareMatrixForGLPK(const Matrix &matrix,
                          LinearProblem &linear_problem);

#endif //CFLEX_FIND_INITIAL_FLUXES_H

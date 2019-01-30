#ifndef CFLEX_METABOLIC_FLUX_ANALYSIS_H
#define CFLEX_METABOLIC_FLUX_ANALYSIS_H

#include "../utilities/reaction_struct.h"
#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/objective_parameters.h"
#include "../math/math_utilites.h"

#include "../alglib/optimization.h"

#include <vector>
#include <string>
#include <random>

std::vector<Flux> EstimateFluxes(ObjectiveParameters *parameters,
                                 const std::vector<FluxVariability> &flux_ranges,
                                 const Matrix &stoichiometry_matrix,
                                 const std::vector<Reaction> &reactions,
                                 const int iteration_total);

void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr);


double SimulateAndGetSSR(alglib::real_1d_array &free_fluxes,
                         alglib::real_1d_array lower_bounds,
                         alglib::real_1d_array upper_bounds);


void FillBoundVectors(alglib::real_1d_array &lower_bounds,
                      alglib::real_1d_array &upper_bounds,
                      const std::vector<FluxVariability> &flux_ranges,
                      const std::vector<Reaction> &reactions,
                      const int nullity);


void GenerateInitialPoints(alglib::real_1d_array &free_fluxes,
                           const alglib::real_1d_array &lower_bounds,
                           const alglib::real_1d_array &upper_bounds,
                           const std::vector<Reaction> &reactions,
                           const int nullity,
                           std::mt19937 &random_source);


int GetMeasurementsCount(ObjectiveParameters *parameters);


std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes,
                                             const Matrix &nullspace,
                                             const std::vector<Reaction> &reactions);

void Fillf0Array(alglib::real_1d_array &residuals,
                 const std::vector<EMUandMID> &simulated_mids,
                 const ObjectiveParameters &parameters);


double GetSSR(const alglib::real_1d_array &residuals, int measurements_count);

#endif //CFLEX_METABOLIC_FLUX_ANALYSIS_H

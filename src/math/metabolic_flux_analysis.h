#ifndef CFLEX_METABOLIC_FLUX_ANALYSIS_H
#define CFLEX_METABOLIC_FLUX_ANALYSIS_H

#include "reaction_struct.h"
#include "EMU.h"
#include "MID.h"
#include "objective_parameters.h"
#include "math_utilites.h"

#include "optimization.h"

#include <vector>
#include <string>
#include <random>
#include <tuple>

std::vector<Flux> EstimateFluxes(ObjectiveParameters &parameters,
                                 const int iteration_total);

int GetMeasurementsCount(const std::vector<Measurement> &measurements);

void FillBoundVectors(alglib::real_1d_array &lower_bounds,
                      alglib::real_1d_array &upper_bounds,
                      const std::vector<Reaction> &reactions,
                      const int nullity);


void GenerateInitialPoints(alglib::real_1d_array &free_fluxes,
                           const alglib::real_1d_array &lower_bounds,
                           const alglib::real_1d_array &upper_bounds,
                           const std::vector<Reaction> &reactions,
                           const int nullity,
                           std::mt19937 &random_source);


void SetOptimizationParameters(alglib::real_1d_array &free_fluxes,
                               alglib::real_1d_array &lower_bounds,
                               alglib::real_1d_array &upper_bounds,
                               int nullity,
                               int measurements_count,
                               alglib::minlmstate &state,
                               alglib::minlmreport &report);


std::tuple<double, std::vector<Flux>> RunOptimization( int measurements_count,
                                                       ObjectiveParameters *parameters,
                                                       alglib::minlmstate &state,
                                                       alglib::minlmreport &report);

void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr);

std::vector<Flux> CalculateAllFluxesFromFree(const alglib::real_1d_array &free_fluxes,
                                             const Matrix &nullspace,
                                             const std::vector<Reaction> &reactions);


Eigen::VectorXd GetEigenVectorFromAlgLibVector(const alglib::real_1d_array &alglib_vector);

void Fillf0Array(alglib::real_1d_array &residuals,
                 const std::vector<EMUandMID> &simulated_mids,
                 const std::vector<Measurement> &measurements);


double GetSSR(const alglib::real_1d_array &residuals, int measurements_count);


#endif //CFLEX_METABOLIC_FLUX_ANALYSIS_H

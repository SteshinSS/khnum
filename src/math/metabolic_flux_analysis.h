#ifndef CFLEX_METABOLIC_FLUX_ANALYSIS_H
#define CFLEX_METABOLIC_FLUX_ANALYSIS_H

#include "../utilities/reaction_struct.h"
#include "../utilities/EMU.h"
#include "../utilities/MID.h"

#include "../alglib/optimization.h"

#include <vector>
#include <string>
#include <map>

std::map<std::string, Flux> EstimateFluxes(const std::vector<EMUNetwork> &networks,
                                           const std::map<std::string, Flux> &initial_fluxes,
                                           const std::map<std::string, FluxVariability> &flux_ranges);


void CalculateResidual(const alglib::real_1d_array &free_fluxes,
                       alglib::real_1d_array &residuals,
                       void *ptr);

#endif //CFLEX_METABOLIC_FLUX_ANALYSIS_H

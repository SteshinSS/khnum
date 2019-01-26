#include "metabolic_flux_analysis.h"

#include "../utilities/reaction_struct.h"
#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../math/calculate_mids.h"

#include <vector>
#include <string>
#include <map>

std::map<std::string, Flux> EstimateFluxes(const std::vector<EMUNetwork> &networks,
                                           const std::map<std::string, Flux> &initial_fluxes,
                                           const std::map<std::string, FluxVariability> &flux_ranges) {

}

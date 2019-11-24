#pragma once

#include "matrix.h"
#include "input_substrate.h"
#include "emu_and_mid.h"
#include "emu.h"
#include "measurement.h"
#include "reaction.h"


namespace khnum {
struct ReactionsName {
    int id;
    std::string name;
};

struct GeneratorParameters {
    std::vector<EmuNetwork> networks;
    std::vector<EmuAndMid> input_mids;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<int> free_flux_id_to_nullspace_position;
    std::vector<int> free_fluxes_id;
};

struct Problem {
    std::vector<ReactionsName> reactions;
    size_t reactions_total;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<Measurement> measurements;
    int measurements_count;
    bool use_analytic_jacobian = false;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    GeneratorParameters simulator_parameters_;
};
} //namesapce khnum
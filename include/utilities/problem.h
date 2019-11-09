#pragma once

#include "matrix.h"
#include "input_substrate.h"
#include "emu_and_mid.h"
#include "emu.h"
#include "measurement.h"
#include "reaction.h"


namespace khnum {

struct SimulatorParameters {
    std::vector<EmuNetwork> networks;
    std::vector<EmuAndMid> input_mids;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<int> id_to_pos;
    std::vector<int> free_fluxes_id;
    bool use_analytic_jacobian = false;
};

struct Problem {
    std::vector<Reaction> reactions;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<Measurement> measurements;
    int measurements_count;
    SimulatorParameters simulator_parameters_;
};
} //namesapce khnum
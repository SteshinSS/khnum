#pragma once

#include "matrix.h"
#include "input_substrate.h"
#include "emu_and_mid.h"
#include "emu.h"
#include "measurement.h"
#include "reaction.h"


namespace khnum {
struct Problem {
    std::vector<Reaction> reactions;
    std::vector<Emu> measured_isotopes;
    Matrix nullspace;
    std::vector<EmuNetwork> networks;
    std::vector<EmuAndMid> input_substrate_mids;
    std::vector<Measurement> measurements;
    int measurements_count;
};
} //namesapce khnum
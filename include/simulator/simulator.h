#pragma once

#include <vector>

#include "utilities/emu_and_mid.h"
#include "utilities/matrix.h"
#include "utilities/reaction.h"
#include "simulator/flux_combination.h"
#include "simulator/simulation_data.h"


namespace khnum {



class Simulator {
public:
    Simulator(const std::vector<SimulatorNetworkData>& networks,
              const std::vector<EmuAndMid>& input_mids,
              const size_t total_mids_to_simulate);

    SimulatorResult CalculateMids(const std::vector<Flux> &fluxes, bool calculate_jacobian);

private:
    const size_t total_networks_;
    const size_t total_free_fluxes_;
    const size_t total_mids_to_simulate_;
    const std::vector<EmuAndMid> input_mids_;
    std::vector<SimulatorNetworkData> networks_;
};
} // namespace khnum
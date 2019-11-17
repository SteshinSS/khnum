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
              const size_t total_free_fluxes);

    SimulatorResult CalculateMids(const std::vector<Flux> &fluxes, bool calculate_jacobian);

private:
    const size_t total_networks_;
    const size_t total_free_fluxes_;
    const std::vector<EmuAndMid> input_mids_;
    std::vector<SimulatorNetworkData> networks_;

    std::vector<EmuAndMid> simulated_mids_;
    std::vector<std::vector<EmuAndMid>> diff_results_;

};
} // namespace khnum
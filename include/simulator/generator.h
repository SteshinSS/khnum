#pragma once

#include "utilities/problem.h"
#include "simulator/new_simulator.h"

namespace khnum {
class SimulatorGenerator {
public:
    SimulatorGenerator(const SimulatorParameters& parameters);

    NewSimulator Generate();

private:
    std::vector<NetworkEmu> InitializeInputEmus(const std::vector<EmuAndMid>& input_mids) const;
    std::vector<Emu> measured_isotopes_;
    std::vector<SimulatorNetworkData> simulator_network_data_;
    std::vector<EmuAndMid> input_mids_;
    int total_free_fluxes_;
};
}
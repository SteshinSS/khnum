#pragma once

#include "utilities/problem.h"
#include "simulator/simulator.h"

namespace khnum {
class SimulatorGenerator {
public:
    SimulatorGenerator(const GeneratorParameters& parameters);

    Simulator Generate();

private:
    std::vector<NetworkEmu> InitializeInputEmus(const std::vector<EmuAndMid>& input_mids) const;
    SimulatorNetworkData FillSimulatorNetworkData(const GeneratorNetworkData& network_data, int network_size) const;

    GeneratorParameters parameters_;

    std::vector<SimulatorNetworkData> simulator_network_data_;
};
}
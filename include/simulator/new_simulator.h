#pragma once

#include <vector>
#include <alglib/ap.h>

#include "utilities/emu_and_mid.h"
#include "utilities/matrix.h"
#include "utilities/reaction.h"
#include "simulator/flux_combination.h"


namespace khnum {
struct SimulatorResult {
    std::vector<EmuAndMid> simulated_mids;
    std::vector<std::vector<EmuAndMid>> diff_results;
};

struct DerivativeData {
    Matrix dA;
    Matrix dB;
    Matrix dY;
};

struct SimulatorNetworkData {
    std::vector<FluxCombination> symbolic_A;
    std::vector<FluxCombination> symbolic_B;
    std::vector<PositionOfSavedEmu> Y_data;
    std::vector<Convolution> convolutions;
    std::vector<int> usefull_emus;
    std::vector<FinalEmu> final_emus;
    Matrix A;
    Matrix B;
    Matrix Y;
    std::vector<DerivativeData> derivatives;
};


class NewSimulator {
public:
    SimulatorResult CalculateMids(const std::vector<Flux> &fluxes);

private:
    const size_t total_networks;
    const size_t total_free_fluxes;
    const size_t total_measurements;
    const std::vector<EmuAndMid> input_mids_;
    std::vector<SimulatorNetworkData> networks_;
};
} // namespace khnum
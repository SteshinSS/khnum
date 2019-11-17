#pragma once

#include <vector>

#include "utilities/matrix.h"
#include "utilities/emu_and_mid.h"
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

struct GeneratorNetworkData {
    int network_num;
    std::vector<Emu> unknown_emus;
    std::vector<Emu> known_emus;
    std::vector<Convolution> convolutions;
    std::vector<FinalEmu> final_emus;
    std::vector<PositionOfSavedEmu> Y_data;
    std::vector<FluxCombination> symbolic_A;
    std::vector<FluxCombination> symbolic_B;
    std::unordered_map<int, int> reaction_to_convolution;
};
}
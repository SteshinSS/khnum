#pragma once

#include <vector>

#include "utilities/matrix.h"
#include "utilities/emu_and_mid.h"
#include "simulator/flux_combination.h"

namespace khnum {

// contains result of simulation
struct SimulatorResult {
    std::vector<EmuAndMid> simulated_mids;

    // contains derivative of i'th mid by v free flux at [v][i]
    std::vector<std::vector<EmuAndMid>> diff_results;
};

struct DerivativeData {
    std::vector<Triplet> symbolic_dA;
    std::vector<Triplet> symbolic_dB;
};

enum class NetworkSize {small, big};

struct SimulatorNetworkData {
    NetworkSize size;
    std::vector<FluxCombination> symbolic_A;
    std::vector<FluxCombination> symbolic_B;
    std::vector<PositionOfSavedEmu> Y_data;
    std::vector<Convolution> convolutions;
    std::vector<int> usefull_emus;
    std::vector<FinalEmu> final_emus;
    size_t A_rows;
    size_t A_cols;
    size_t B_rows;
    size_t B_cols;
    size_t Y_rows;
    size_t Y_cols;
    std::vector<DerivativeData> derivatives;
};

struct GeneratorNetworkData {
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
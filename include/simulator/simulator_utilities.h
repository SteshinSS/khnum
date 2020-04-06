#pragma once

#include <vector>

#include "utilities/matrix.h"
#include "simulator/flux_combination.h"

namespace khnum {
namespace simulator_utilities {
void FillSmallFluxMatrix(const std::vector<FluxCombination> &symbolic_matrix,
                         const std::vector<Flux> &fluxes,
                         Matrix &matrix_out);

void FillBigFluxMatrix(const std::vector<FluxCombination> &symbolic_matrix,
                       const std::vector<Flux> &fluxes,
                       std::vector<Triplet> &triplet_out);

void FillYMatrix(const std::vector<PositionOfSavedEmu> &Y_data,
                 const std::vector<EmuAndMid> &input_mids,
                 const std::vector<std::vector<Mid>> &saved_mids,
                 const std::vector<Convolution> &convolutions,
                 Matrix &Y_out);

void SaveNewEmus(const Matrix& X,
                 const std::vector<int>& usefull_emus,
                 const std::vector<FinalEmu>& final_emus,
                 std::vector<Mid>& saved_mids_out,
                 std::vector<EmuAndMid>& result_out,
                 std::vector<double>& sums_out);

void SaveNewDiffEmus(const Matrix& X,
                     const std::vector<int>& usefull_emus,
                     const std::vector<FinalEmu>& final_emus,
                     const std::vector<EmuAndMid> result,
                     const std::vector<double>& sums,
                     std::vector<Mid>& saved_mids_out,
                     std::vector<EmuAndMid>& diff_result_out);

Matrix GetCorrectedDiffMid(const Matrix& correction_matrix,
                           Matrix mid,
                           const Matrix& diff_mid,
                           double sum);

Mid ConvolvePartialDiff(const Convolution& convolution,
                        const std::vector<std::vector<Mid>>& known_d_mids,
                        const std::vector<EmuAndMid>& input_mids,
                        const std::vector<std::vector<Mid>>& saved_mids,
                        size_t mid_size,
                        size_t diff_position);

void FillDiffYMatrix(const std::vector<PositionOfSavedEmu>& Y_data,
                     const std::vector<std::vector<Mid>>& known_d_mids,
                     const std::vector<Convolution>& convolutions,
                     const std::vector<EmuAndMid>& input_mids,
                     const std::vector<std::vector<Mid>>& saved_mids,
                     Matrix& Y_out);

}
}
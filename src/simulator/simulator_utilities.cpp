#include "simulator/simulator_utilities.h"

namespace khnum {
namespace simulator_utilities {
void FillFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                    const std::vector<Flux>& fluxes,
                    Matrix& matrix_out) {
    for (const FluxCombination& combination : symbolic_matrix) {
        double value = 0.0;
        for (const FluxAndCoefficient& flux : combination.fluxes) {
            value += flux.coefficient * fluxes[flux.id];
        }
        matrix_out(combination.i, combination.j) = value;
    }
}

void FillYMatrix(const std::vector<PositionOfSavedEmu>& Y_data,
                 const std::vector<EmuAndMid>& input_mids,
                 const std::vector<std::vector<Mid>>& saved_mids,
                 const std::vector<Convolution>& convolutions,
                 Matrix& Y_out) {
    for (size_t i = 0; i < Y_data.size(); ++i) {
        const PositionOfSavedEmu& known_emu = Y_data[i];
        Mid mid;
        if (known_emu.network == -1) {
            mid = input_mids[known_emu.position].mid;
        } else {
            mid = saved_mids[known_emu.network][known_emu.position];
        }
        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y_out(i, mass_shift) = mid[mass_shift];
        }
    }

    for (size_t i = 0; i < convolutions.size(); ++i) {
        const Convolution& convolution = convolutions[i];
        Mid mid(1, 1.0); // MID = [1.0]
        for (const PositionOfSavedEmu& emu : convolution.elements) {
            if (emu.network == -1) {
                mid = mid * input_mids[emu.position].mid;
            } else {
                mid = mid * saved_mids[emu.network][emu.position];
            }
        }
        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y_out(i + Y_data.size(), mass_shift) = mid[mass_shift];
        }
    }
}

void SaveNewEmus(const Matrix& X,
                 const std::vector<int>& usefull_emus,
                 const std::vector<FinalEmu>& final_emus,
                 std::vector<Mid>& saved_mids_out,
                 std::vector<EmuAndMid>& result_out) {
    for (int position : usefull_emus) {
        Mid new_mid;
        for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
            new_mid.push_back(X(position, mass_shift));
        }
        saved_mids_out.push_back(new_mid);
    }

    for (const FinalEmu& final_emu : final_emus) {
        Mid result_mid;
        for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
            result_mid.push_back(X(final_emu.order_in_X, mass_shift));
        }
        EmuAndMid result_emu;
        result_emu.emu = final_emu.emu;
        result_emu.mid = result_mid;
        result_out[final_emu.position_in_result] = result_emu;
    }
}

Mid ConvolvePartialDiff(const Convolution& convolution,
                        const std::vector<std::vector<Mid>>& known_d_mids,
                        const std::vector<EmuAndMid>& input_mids,
                        const std::vector<std::vector<Mid>>& saved_mids,
                        int mid_size,
                        int diff_position) {
    Mid mid_part(1, 1.0);
    for (int j = 0; j < convolution.elements.size(); ++j) {
        const PositionOfSavedEmu& emu = convolution.elements[j];
        if (diff_position == j) {
            if (emu.network == -1) {
                return std::vector<double> (mid_size, 0.0);
            }
            mid_part = mid_part * known_d_mids[emu.network][emu.position];
        } else {
            if (emu.network == -1) {
                mid_part = mid_part * input_mids[emu.position].mid;
            } else {
                mid_part = mid_part * saved_mids[emu.network][emu.position];
            }

        }
    }
    return mid_part;
}

void FillDiffYMatrix(const std::vector<PositionOfSavedEmu>& Y_data,
                     const std::vector<std::vector<Mid>>& known_d_mids,
                     const std::vector<Convolution>& convolutions,
                     const std::vector<EmuAndMid>& input_mids,
                     const std::vector<std::vector<Mid>>& saved_mids,
                     Matrix& Y_out) {
    for (size_t i = 0; i < Y_data.size(); ++i) {
        const PositionOfSavedEmu known_emu = Y_data[i];
        Mid mid;
        if (known_emu.network == -1) {
            // Do nothing. Because:
            // mid = std::vector<double> (Y_out.cols(), 0.0);
        } else {
            mid = known_d_mids[known_emu.network][known_emu.position];
        }
        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y_out(i, mass_shift) = mid[mass_shift];
        }
    }

    size_t position = Y_data.size();
    for (const Convolution& convolution : convolutions) {
        Mid mid = std::vector<double> (Y_out.cols(), 0.0);
        for (int diff_position = 0; diff_position < convolution.elements.size(); ++diff_position) {
            Mid mid_part = ConvolvePartialDiff(convolution, known_d_mids, input_mids, saved_mids, Y_out.cols(), diff_position);
            for (int j = 0; j < mid.size(); ++j) {
                mid[j] += mid_part[j];
            }
        }

        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y_out(position, mass_shift) = mid[mass_shift];
        }
        ++position;
    }
}
}
}
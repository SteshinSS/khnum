#include "simulator/simulator_utilities.h"

namespace khnum {
namespace simulator_utilities {
void FillSmallFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
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

void FillBigFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                       const std::vector<Flux>& fluxes,
                       std::vector<Triplet> &triplets_out) {
    for (const FluxCombination& combination : symbolic_matrix) {
        double value = 0.0;
        for (const FluxAndCoefficient& flux : combination.fluxes) {
            value += flux.coefficient * fluxes[flux.id];
        }
        triplets_out.push_back(Triplet(combination.i, combination.j, value));
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
                 std::vector<EmuAndMid>& result_out,
                 std::vector<double>& sums_out) {
    for (int position : usefull_emus) {
        Mid new_mid;
        for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
            new_mid.push_back(X(position, mass_shift));
        }
        saved_mids_out.push_back(new_mid);
    }

    for (const FinalEmu& final_emu : final_emus) {
        Mid result_mid;
        if (final_emu.correction_matrix.rows() > 0) {
            Matrix original_mid = X.row(final_emu.order_in_X);
            Matrix corrected_mid = final_emu.correction_matrix * original_mid.transpose();
            double sum = 0.0;
            for (int i = 0; i < corrected_mid.rows(); ++i) {
                sum += corrected_mid(i, 0);
            }
            sums_out[final_emu.position_in_result] = sum;
            for (int mass_shift = 0; mass_shift < corrected_mid.rows(); ++mass_shift) {
                result_mid.push_back(corrected_mid(mass_shift, 0) / sum);
            }
        } else {
            for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
                result_mid.push_back(X(final_emu.order_in_X, mass_shift));
            }
        }

        EmuAndMid result_emu;
        result_emu.emu = final_emu.emu;
        result_emu.mid = result_mid;
        result_out[final_emu.position_in_result] = result_emu;
    }
}

void SaveNewDiffEmus(const Matrix& X,
                     const std::vector<int>& usefull_emus,
                     const std::vector<FinalEmu>& final_emus,
                     const std::vector<EmuAndMid> result,
                     const std::vector<double>& sums,
                     std::vector<Mid>& saved_mids_out,
                     std::vector<EmuAndMid>& diff_result_out) {
    for (int position : usefull_emus) {
        Mid new_mid;
        for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
            new_mid.push_back(X(position, mass_shift));
        }
        saved_mids_out.push_back(new_mid);
    }

    for (const FinalEmu& final_emu : final_emus) {
        Mid result_mid;
        if (final_emu.correction_matrix.rows() > 0) {
            Matrix mid(result[final_emu.position_in_result].mid.size(), 1);
            for (int i = 0; i < result[final_emu.position_in_result].mid.size(); ++i) {
                mid(i, 0) = result[final_emu.position_in_result].mid[i];
            }

            Matrix diff_mid = X.row(final_emu.order_in_X);
            Matrix corrected_mid = GetCorrectedDiffMid(final_emu.correction_matrix,
                                                       mid,
                                                       diff_mid.transpose(),
                                                       sums[final_emu.position_in_result]);


            for (int mass_shift = 0; mass_shift < corrected_mid.rows(); ++mass_shift) {
                result_mid.push_back(corrected_mid(mass_shift, 0));
            }
        } else {
            for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
                result_mid.push_back(X(final_emu.order_in_X, mass_shift));
            }
        }

        EmuAndMid result_emu;
        result_emu.emu = final_emu.emu;
        result_emu.mid = result_mid;
        diff_result_out[final_emu.position_in_result] = result_emu;
    }
}

Matrix GetCorrectedDiffMid(const Matrix& correction_matrix,
                           Matrix mid,
                           const Matrix& diff_mid,
                           double sum) {
    Matrix corrected_diff_mid = correction_matrix * diff_mid;
    double diff_sum = 0.0;
    for (int i = 0; i < corrected_diff_mid.rows(); ++i) {
        diff_sum += corrected_diff_mid(i, 0);
    }
    mid *= diff_sum;
    corrected_diff_mid -= mid;
    corrected_diff_mid /= sum;
    return corrected_diff_mid;
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
            // and Y_out(i, ...) is already zero
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
        for (size_t diff_position = 0; diff_position < convolution.elements.size(); ++diff_position) {
            Mid mid_part = ConvolvePartialDiff(convolution,
                                               known_d_mids,
                                               input_mids,
                                               saved_mids,
                                               Y_out.cols(),
                                               diff_position);
            for (size_t j = 0; j < mid.size(); ++j) {
                mid[j] += mid_part[j];
            }
        }

        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y_out(position, mass_shift) = mid[mass_shift];
        }
        ++position;
    }
}

Mid ConvolvePartialDiff(const Convolution& convolution,
                        const std::vector<std::vector<Mid>>& known_d_mids,
                        const std::vector<EmuAndMid>& input_mids,
                        const std::vector<std::vector<Mid>>& saved_mids,
                        size_t mid_size,
                        size_t diff_position) {
    Mid mid_part(1, 1.0);
    for (size_t i = 0; i < convolution.elements.size(); ++i) {
        const PositionOfSavedEmu& emu = convolution.elements[i];
        if (diff_position == i) {
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
} // namespace simulator_utilities
} // namespace khnum
#include "simulator/simulator.h"

#include <vector>
#include <iostream>

#include "simulator/simulation_data.h"
#include "simulator/simulator_utilities.h"

namespace khnum {


SimulatorResult Simulator::CalculateMids(const std::vector<Flux> &fluxes, bool calculate_jacobian) {
    SimulatorResult result;
    std::vector<double> sums(total_mids_to_simulate_);
    std::vector<EmuAndMid>& simulated_mids = result.simulated_mids;
    simulated_mids.resize(total_mids_to_simulate_);

    std::vector<std::vector<EmuAndMid>>& diff_results = result.diff_results;
    diff_results.resize(total_free_fluxes_);
    for (std::vector<EmuAndMid>& vec : diff_results) {
        vec.resize(total_mids_to_simulate_);
    }

    std::vector<std::vector<Mid>> saved_mids(total_networks_);
    // contains MID's derivative at [free_flux][network][i]
    std::vector<std::vector<std::vector<Mid>>> saved_diff_mids(total_free_fluxes_);
    for (auto& vec : saved_diff_mids) {
        vec.resize(total_networks_);
    }
    for (size_t network_num = 0; network_num < total_networks_; ++network_num) {
        const SimulatorNetworkData &network = networks_[network_num];
        if (network.size == NetworkSize::small) {
            Matrix A = Matrix::Zero(network.A_rows, network.A_cols);
            simulator_utilities::FillSmallFluxMatrix(network.symbolic_A, fluxes, A);

            Matrix B = Matrix::Zero(network.B_rows, network.B_cols);
            simulator_utilities::FillSmallFluxMatrix(network.symbolic_B, fluxes, B);

            Matrix Y = Matrix::Zero(network.Y_rows, network.Y_cols);
            simulator_utilities::FillYMatrix(network.Y_data,input_mids_, saved_mids,
                                             network.convolutions, Y);

            const Matrix BY = B * Y;
            const Matrix X = A.householderQr().solve(BY);
            simulator_utilities::SaveNewEmus(X, network.usefull_emus, network.final_emus, saved_mids[network_num], simulated_mids,
                                             sums);
            if (!calculate_jacobian) {
                continue;
            }
            Matrix dY = Matrix::Zero(network.Y_rows, network.Y_cols);
            Matrix dA = Matrix::Zero(network.A_rows, network.A_cols);
            Matrix dB = Matrix::Zero(network.B_rows, network.B_cols);
            for (size_t flux = 0; flux < total_free_fluxes_; ++flux) {
                const DerivativeData &derivatives = network.derivatives.at(flux);

                dY.setZero();
                simulator_utilities::FillDiffYMatrix(network.Y_data, saved_diff_mids[flux], network.convolutions,
                                                     input_mids_, saved_mids, dY);

                dA.setZero();
                for (const SparseElement &element : derivatives.symbolic_dA) {
                    dA(element.i, element.j) = element.value;
                }

                dB.setZero();
                for (const SparseElement &element : derivatives.symbolic_dB) {
                    dB(element.i, element.j) = element.value;
                }

                // Right Part of A * dX = (...)
                Matrix RightPart = dB * Y + B * dY - dA * X;
                Matrix dX = A.householderQr().solve(RightPart);
                simulator_utilities::SaveNewDiffEmus(dX, network.usefull_emus, network.final_emus, simulated_mids,
                                                     sums, saved_diff_mids[flux][network_num], diff_results[flux]);
            }
        } else {
            SparseMatrix A(network.A_rows, network.A_cols);
            std::vector<Triplet> A_triplets;
            simulator_utilities::FillBigFluxMatrix(network.symbolic_A, fluxes, A_triplets);
            A.setFromTriplets(A_triplets.begin(), A_triplets.end());

            SparseMatrix B(network.B_rows, network.B_cols);
            std::vector<Triplet> B_triplets;
            simulator_utilities::FillBigFluxMatrix(network.symbolic_B, fluxes, B_triplets);
            B.setFromTriplets(B_triplets.begin(), B_triplets.end());

            Matrix Y = Matrix::Zero(network.Y_rows, network.Y_cols);
            simulator_utilities::FillYMatrix(network.Y_data,input_mids_, saved_mids,
                                             network.convolutions, Y);

            const Matrix BY = B * Y;
            Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<SparseMatrix::Scalar>> solver;
            solver.compute(A);
            const Matrix X = solver.solve(BY);
            if (solver.info() != Eigen::Success) {
                std::cout << "NOT SUCCESS " << solver.info() << std::endl;
                std::cout << solver.error() << std::endl;
                std::cout << solver.iterations() << std::endl;
            }

            simulator_utilities::SaveNewEmus(X, network.usefull_emus, network.final_emus, saved_mids[network_num], simulated_mids,
                                             sums);
            if (!calculate_jacobian) {
                continue;
            }
        }


    }
    return result;
}

Simulator::Simulator(const std::vector<SimulatorNetworkData>& networks,
                     const std::vector<EmuAndMid>& input_mids,
                     const size_t total_mids_to_simulate) :
                                            total_networks_{networks.size()},
                                            total_free_fluxes_{networks[0].derivatives.size()},
                                            total_mids_to_simulate_{total_mids_to_simulate},
                                            input_mids_{input_mids},
                                            networks_{networks} {
}
} // namespace khnum
#include "simulator/simulator.h"

#include <vector>
#include <iostream>
#include <chrono>

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
    size_t total_big_networks = 0;
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
            auto A_decomposition = A.householderQr();
            const Matrix X = A_decomposition.solve(BY);
            simulator_utilities::SaveNewEmus(X, network.usefull_emus, network.final_emus, saved_mids[network_num], simulated_mids,
                                             sums);
            if (!calculate_jacobian) {
                continue;
            }
            Matrix dY = Matrix::Zero(network.Y_rows, network.Y_cols);
            for (size_t flux = 0; flux < total_free_fluxes_; ++flux) {
                const DerivativeData &derivatives = network.derivatives.at(flux);

                dY.setZero();
                simulator_utilities::FillDiffYMatrix(network.Y_data, saved_diff_mids[flux], network.convolutions,
                                                     input_mids_, saved_mids, dY);

                // Right Part of A * dX = (...)
                Matrix RightPart = derivatives.dB_small * Y + B * dY - derivatives.dA_small * X;
                Matrix dX = A_decomposition.solve(RightPart);
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
            Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<SparseMatrix::Scalar>> &solver = solvers_[total_big_networks++];
            solver.factorize(A);
            const Matrix guess = Matrix::Constant(network.A_rows, network.Y_cols, 1.0 / network.Y_cols);

            const Matrix X = solver.solveWithGuess(BY, guess);
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
            Matrix dY = Matrix::Zero(network.Y_rows, network.Y_cols);
            for (size_t flux = 0; flux < total_free_fluxes_; ++flux) {
                const DerivativeData &derivatives = network.derivatives.at(flux);

                dY.setZero();
                simulator_utilities::FillDiffYMatrix(network.Y_data, saved_diff_mids[flux], network.convolutions,
                                                     input_mids_, saved_mids, dY);

                // Right Part of A * dX = (...)
                Matrix RightPart = derivatives.dB_big * Y - derivatives.dA_big * X + B * dY;
                Matrix dX = solver.solve(RightPart);


                if (solver.info() != Eigen::Success) {
                    std::cout << "NOT SUCCESS " << solver.info() << std::endl;
                    std::cout << solver.error() << std::endl;
                    std::cout << solver.iterations() << std::endl;
                }
                simulator_utilities::SaveNewDiffEmus(dX, network.usefull_emus, network.final_emus, simulated_mids,
                                                     sums, saved_diff_mids[flux][network_num], diff_results[flux]);
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
                                            networks_{networks},
                                            solvers_{networks_.size()}{
    int total_big_networks = 0;
    for (const SimulatorNetworkData &network : networks) {
        if (network.size == NetworkSize::big) {
            SparseMatrix A(network.A_rows, network.A_cols);
            std::vector<Triplet> A_triplets;
            for (const FluxCombination &fc : network.symbolic_A) {
                A_triplets.push_back(Triplet(fc.i, fc.j, 1));
            }
            A.setFromTriplets(A_triplets.begin(), A_triplets.end());
            solvers_[total_big_networks].analyzePattern(A);
            ++total_big_networks;
        }
    }
}
} // namespace khnum
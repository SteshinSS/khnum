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
        SimulatorNetworkData& network = networks_[network_num];

        simulator_utilities::FillFluxMatrix(network.symbolic_A, fluxes, network.A);
        simulator_utilities::FillFluxMatrix(network.symbolic_B, fluxes, network.B);
        simulator_utilities::FillYMatrix(network.Y_data,input_mids_, saved_mids,
                                         network.convolutions, network.Y);

        const Matrix BY = network.B * network.Y;
        const Matrix X = network.A.householderQr().solve(BY);
        simulator_utilities::SaveNewEmus(X, network.usefull_emus, network.final_emus, saved_mids[network_num], simulated_mids,
                                        sums);

        if (!calculate_jacobian) {
            continue;
        }
        for (size_t flux = 0; flux < total_free_fluxes_; ++flux) {
            DerivativeData& derivatives = network.derivatives[flux];
            simulator_utilities::FillDiffYMatrix(network.Y_data, saved_diff_mids[flux], network.convolutions,
                                                 input_mids_, saved_mids, derivatives.dY);
            // Right Part of A * dX = (...)
            Matrix RightPart = derivatives.dB * network.Y + network.B * derivatives.dY - derivatives.dA * X;
            Matrix dX = network.A.householderQr().solve(RightPart);
            simulator_utilities::SaveNewDiffEmus(dX, network.usefull_emus, network.final_emus, simulated_mids,
                                             sums, saved_diff_mids[flux][network_num], diff_results[flux]);
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
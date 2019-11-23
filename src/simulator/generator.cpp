#include "simulator/generator.h"
#include "simulator/generator_utilites.h"
#include "simulator/simulator.h"
#include "simulator/simulation_data.h"

namespace khnum {
SimulatorGenerator::SimulatorGenerator(const GeneratorParameters& parameters) {
    parameters_ = parameters;

    std::vector<std::vector<int>> usefull_emus(parameters.networks.size());
    std::vector<NetworkEmu> all_known_emus = InitializeInputEmus(parameters.input_mids);
    for (size_t network_num = 0; network_num < parameters.networks.size(); ++network_num) {
        GeneratorNetworkData network_data;
        const std::vector<EmuReaction>& reactions = parameters.networks[network_num];

        generator_utilites::FillEmuLists(reactions, all_known_emus, network_data, usefull_emus);
        generator_utilites::CreateSymbolicMatrices(reactions, network_data);
        generator_utilites::FillFinalEmu(parameters.measured_isotopes, network_data);
        generator_utilites::InsertIntoAllKnownEmus(network_data.unknown_emus, network_num, all_known_emus);

        int network_size = generator_utilites::FindNetworkSize(reactions);
        simulator_network_data_.emplace_back(FillSimulatorNetworkData(network_data, network_size));
    }

    for (size_t network_num = 0; network_num < parameters.networks.size(); ++network_num) {
        for (size_t position : usefull_emus[network_num]) {
            simulator_network_data_[network_num].usefull_emus.push_back(position);
        }
    }
}

Simulator SimulatorGenerator::Generate() {
    return Simulator(simulator_network_data_, parameters_.input_mids, parameters_.measured_isotopes.size());
}

std::vector<NetworkEmu> SimulatorGenerator::InitializeInputEmus(const std::vector<EmuAndMid>& input_mids) const {
    std::vector<NetworkEmu> all_known_emus;
    for (size_t i = 0; i < input_mids.size(); ++i) {
        NetworkEmu input_emu;
        input_emu.emu = input_mids[i].emu;
        input_emu.network = -1;
        input_emu.order_in_usefull_emus = i;
        input_emu.order_in_X = i;
        input_emu.is_usefull = true;
        all_known_emus.push_back(input_emu);
    }

    return all_known_emus;
}

SimulatorNetworkData SimulatorGenerator::FillSimulatorNetworkData(const GeneratorNetworkData& network_data, int network_size) const {
    SimulatorNetworkData simulator_network_data;
    simulator_network_data.symbolic_A = network_data.symbolic_A;
    simulator_network_data.symbolic_B = network_data.symbolic_B;
    simulator_network_data.Y_data = network_data.Y_data;
    simulator_network_data.convolutions = network_data.convolutions;

    simulator_network_data.final_emus = network_data.final_emus;
    simulator_network_data.A = Matrix::Zero(network_data.unknown_emus.size(), network_data.unknown_emus.size());
    simulator_network_data.B = Matrix::Zero(network_data.unknown_emus.size(), network_data.known_emus.size() + network_data.convolutions.size());

    simulator_network_data.Y = Matrix::Zero(network_data.known_emus.size() + network_data.convolutions.size(), network_size + 1);

    std::vector<DerivativeData> derivatives(parameters_.free_fluxes_id.size());
    size_t position = 0;
    for (int id : parameters_.free_fluxes_id) {
        derivatives[position].dA = generator_utilites::GenerateDiffFluxMatrix(network_data.symbolic_A,
                                                                              network_data.unknown_emus.size(),
                                                                              network_data.unknown_emus.size(),
                                                                              id, position, parameters_.free_flux_id_to_nullspace_position, parameters_.nullspace);

        derivatives[position].dB = generator_utilites::GenerateDiffFluxMatrix(network_data.symbolic_B,
                                                                              network_data.unknown_emus.size(),
                                                                              network_data.known_emus.size() + network_data.convolutions.size(),
                                                                              id, position, parameters_.free_flux_id_to_nullspace_position, parameters_.nullspace);

        derivatives[position].dY = Matrix::Zero(network_data.known_emus.size() + network_data.convolutions.size(), network_size + 1);
        ++position;
    }

    simulator_network_data.derivatives = derivatives;
    return simulator_network_data;
}



}

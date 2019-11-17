#include "simulator/generator.h"
#include "simulator/generator_utilites.h"
#include "simulator/new_simulator.h"

namespace khnum {
SimulatorGenerator::SimulatorGenerator(const SimulatorParameters& parameters) {
    std::vector<SimulatorNetworkData> networks(parameters.networks.size());
    total_free_fluxes_ = parameters.free_fluxes_id.size();
    input_mids_ = parameters.input_mids;
    measured_isotopes_ = parameters.measured_isotopes;

    std::vector<std::vector<int>> usefull_emus(parameters.networks.size());
    std::vector<NetworkEmu> all_known_emus = InitializeInputEmus(input_mids_);
    for (int network = 0; network < parameters.networks.size(); ++network) {
        generator_utilites::GeneratorNetworkData network_data;
        network_data.network_num = network;
        const std::vector<EmuReaction>& reactions = parameters.networks[network];
        generator_utilites::FillEmuLists(reactions, all_known_emus, network_data, usefull_emus);
        generator_utilites::CreateSymbolicMatrices(reactions, network_data);
        generator_utilites::FillFinalEmu(measured_isotopes_, network_data);
        generator_utilites::InsertIntoAllKnownEmus(network_data.unknown_emus, network, all_known_emus);

        SimulatorNetworkData simulator_network_data;
        simulator_network_data.symbolic_A = network_data.symbolic_A;
        simulator_network_data.symbolic_B = network_data.symbolic_B;
        simulator_network_data.Y_data = network_data.Y_data;
        simulator_network_data.convolutions = network_data.convolutions;

        simulator_network_data.final_emus = network_data.final_emus;
        simulator_network_data.A = Matrix::Zero(network_data.unknown_emus.size(), network_data.unknown_emus.size());
        simulator_network_data.B = Matrix::Zero(network_data.unknown_emus.size(), network_data.known_emus.size() + network_data.convolutions.size());

        int network_size = generator_utilites::FindNetworkSize(reactions);
        simulator_network_data.Y = Matrix::Zero(network_data.known_emus.size() + network_data.convolutions.size(), network_size + 1);

        std::vector<DerivativeData> derivatives(parameters.free_fluxes_id.size());
        int position = 0;
        for (int id : parameters.free_fluxes_id) {
            derivatives[position].dA = generator_utilites::GenerateDiffFluxMatrix(network_data.symbolic_A,
                                               network_data.unknown_emus.size(),
                                               network_data.unknown_emus.size(),
                                               id, position, parameters.id_to_pos, parameters.nullspace);

            derivatives[position].dB = generator_utilites::GenerateDiffFluxMatrix(network_data.symbolic_B,
                                                                   network_data.unknown_emus.size(),
                                                                   network_data.known_emus.size() + network_data.convolutions.size(),
                                                                   id, position, parameters.id_to_pos, parameters.nullspace);

            derivatives[position].dY = Matrix::Zero(network_data.known_emus.size() + network_data.convolutions.size(), network_size + 1);
            ++position;
        }

        simulator_network_data.derivatives = derivatives;
        simulator_network_data_.push_back(simulator_network_data);
    }

    for (int network = 0; network < parameters.networks.size(); ++network) {
        for (int position : usefull_emus[network]) {
            simulator_network_data_[network].usefull_emus.push_back(position);
        }
    }
}



NewSimulator SimulatorGenerator::Generate() {
    return NewSimulator(simulator_network_data_, input_mids_, total_free_fluxes_);
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



}

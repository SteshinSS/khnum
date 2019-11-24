#include "simulator/generator_utilites.h"

#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>

#include "utilities/emu.h"
#include "simulator/flux_combination.h"
#include "simulator/simulation_data.h"
#include "utilities/matrix.h"


namespace khnum {
namespace generator_utilites {

void FillEmuLists(const std::vector<EmuReaction>& reactions,
                  std::vector<NetworkEmu>& all_known_emus,
                  GeneratorNetworkData& network_data,
                  std::vector<std::vector<int>>& usefull_emus) {
    std::set<Emu> seen_emus;
    size_t reaction_num = 0;
    for (const EmuReaction &reaction : reactions) {
        if (reaction.left.size() == 1) {
            if (seen_emus.find(reaction.left[0].emu) == seen_emus.end()) {
                CheckAndInsertEmu(reaction.left[0].emu, all_known_emus, network_data, usefull_emus);
                seen_emus.insert(reaction.left[0].emu);
            }
        } else {
            Convolution convolution = ConvolveReaction(reaction, all_known_emus, usefull_emus);

            std::vector<Convolution>& convolutions = network_data.convolutions;
            auto already_has = std::find(convolutions.begin(), convolutions.end(), convolution);
            if (already_has == convolutions.end()) {
                convolutions.push_back(convolution);
            } else {
                throw std::runtime_error("wtf");
            }
            network_data.reaction_to_convolution[reaction_num] = convolutions.size() - 1;
        }

        if (seen_emus.find(reaction.right.emu) == seen_emus.end()) {
            network_data.unknown_emus.push_back(reaction.right.emu);
            seen_emus.insert(reaction.right.emu);
        }
        ++reaction_num;
    }
}


void CheckAndInsertEmu(const Emu &emu,
                       std::vector<NetworkEmu>& all_known_emus,
                       GeneratorNetworkData& network_data,
                       std::vector<std::vector<int>>& usefull_emus) {
    auto it = std::find_if(all_known_emus.begin(),
                           all_known_emus.end(),
                           [&emu](const NetworkEmu& known_emu) {
                               return known_emu.emu == emu;
                           });

    if (it != all_known_emus.end()) {
        if (!it->is_usefull) {
            it->is_usefull = true;
            usefull_emus[it->network].push_back(it->order_in_X);
            it->order_in_usefull_emus = usefull_emus[it->network].size() - 1;
        }
        network_data.known_emus.push_back(it->emu);
        network_data.Y_data.push_back({it->network, it->order_in_usefull_emus});
    } else {
        network_data.unknown_emus.push_back(emu);
    }
}


Convolution ConvolveReaction(const EmuReaction& reaction,
                             std::vector<NetworkEmu>& all_known_emus,
                             std::vector<std::vector<int>>& usefull_emus) {
    Convolution convolution;
    convolution.flux_id = reaction.id;
    for (const EmuSubstrate& emu : reaction.left) {
        auto it = std::find_if(all_known_emus.begin(),
                               all_known_emus.end(),
                               [&emu](const NetworkEmu& known_emu) {
                                   return known_emu.emu == emu.emu;
                               });

        if (!it->is_usefull) {
            it->is_usefull = true;
            usefull_emus[it->network].push_back(it->order_in_X);
            it->order_in_usefull_emus = usefull_emus[it->network].size() - 1;
        }

        PositionOfSavedEmu new_known_emu;
        new_known_emu.network = it->network;
        new_known_emu.position = it->order_in_usefull_emus;
        convolution.elements.emplace_back(new_known_emu);
    }

    return convolution;
}


void CreateSymbolicMatrices(const std::vector<EmuReaction>& reactions,
                            GeneratorNetworkData& network_data) {
    const std::vector<Emu>& unknown_emus = network_data.unknown_emus;
    const std::vector<Emu> known_emus = network_data.known_emus;
    const std::vector<Convolution> convolutions = network_data.convolutions;

    const int unknown_size = unknown_emus.size();
    const int known_size = known_emus.size() + convolutions.size();

    // A[unknown_size][unknown_size]
    std::vector<std::vector<FluxCombination>> A(unknown_size, std::vector<FluxCombination>(unknown_size));

    // B[unknown_size][known_size]
    std::vector<std::vector<FluxCombination>> B(unknown_size, std::vector<FluxCombination>(known_size));

    size_t reaction_num = 0;
    for (const EmuReaction &reaction : reactions) {
        size_t position_of_product = FindUnknownEmuPosition(reaction.right.emu, unknown_emus);

        {
            FluxAndCoefficient product;
            product.coefficient = -reaction.rate;
            product.id = reaction.id;
            A[position_of_product][position_of_product].fluxes.emplace_back(product);
        }

        if (reaction.left.size() > 1) {
            int position_of_convolution = network_data.reaction_to_convolution.at(reaction_num);
            FluxAndCoefficient convolution_element;
            convolution_element.coefficient = -reaction.rate;
            convolution_element.id = reaction.id;
            B[position_of_product][known_emus.size() + position_of_convolution].fluxes.emplace_back(convolution_element);

        } else {
            EmuSubstrate substrate = reaction.left[0];

            // returns -1 if emu is unknown
            int position_of_substrate = FindKnownEmuPosition(substrate.emu, known_emus);
            if (position_of_substrate == -1) {
                FluxAndCoefficient substrate_element;
                substrate_element.coefficient = reaction.rate;
                substrate_element.id = reaction.id;
                position_of_substrate = FindUnknownEmuPosition(substrate.emu, unknown_emus);
                A[position_of_product][position_of_substrate].fluxes.emplace_back(substrate_element);
            } else {
                FluxAndCoefficient substrate_element;
                substrate_element.coefficient = -reaction.rate;
                substrate_element.id = reaction.id;
                B[position_of_product][position_of_substrate].fluxes.emplace_back(substrate_element);
            }
        }
        ++reaction_num;
    }
    ConvertToSparseMatrix(A, network_data.symbolic_A);
    ConvertToSparseMatrix(B, network_data.symbolic_B);
}



int FindUnknownEmuPosition(const Emu &emu,
                           const std::vector<Emu>& unknown_emus) {
    auto position = find(unknown_emus.begin(),
                         unknown_emus.end(),
                         emu);

    return position - unknown_emus.begin();
}


int FindKnownEmuPosition(const Emu &emu,
                         const std::vector<Emu>& known_emus) {
    auto position = find_if(known_emus.begin(),
                            known_emus.end(),
                            [&emu](const Emu &known_emu) {
                                return known_emu == emu;
                            });

    if (position != known_emus.end()) {
        return position - known_emus.begin();
    } else {
        return -1;
    }
}

void ConvertToSparseMatrix(const std::vector<std::vector<FluxCombination>>& dense_matrix,
                           std::vector<FluxCombination>& sparse_matrix) {
    for (size_t i = 0; i < dense_matrix.size(); ++i) {
        for (size_t j = 0; j < dense_matrix[0].size(); ++j) {
            if (!dense_matrix[i][j].fluxes.empty()) {
                FluxCombination matrix_element = dense_matrix[i][j];
                matrix_element.i = i;
                matrix_element.j = j;
                sparse_matrix.emplace_back(matrix_element);
            }
        }
    }
}


void FillFinalEmu(const std::vector<Measurement>& measured_isotopes,
                  GeneratorNetworkData& network_data) {
    const std::vector<Emu>& unknown_emus = network_data.unknown_emus;
    for (size_t i = 0; i < unknown_emus.size(); ++i) {
        Emu emu = unknown_emus[i];
        auto it = std::find_if(measured_isotopes.begin(),
                               measured_isotopes.end(),
                               [&emu](const Measurement& measured_emu) {
                                   return measured_emu.emu == emu;
                               });

        if (it != measured_isotopes.end()) {
            FinalEmu final_emu;
            final_emu.emu = emu;
            final_emu.order_in_X = i;
            final_emu.position_in_result = it - measured_isotopes.begin();
            final_emu.correction_matrix = it->correction_matrix;
            network_data.final_emus.push_back(final_emu);
        }
    }
}

void InsertIntoAllKnownEmus(const std::vector<Emu>& unknown_emus,
                            int network_num,
                            std::vector<NetworkEmu>& all_known_emus) {
    for (size_t i = 0; i < unknown_emus.size(); ++i) {
        NetworkEmu new_emu;
        new_emu.order_in_X = i;
        new_emu.emu = unknown_emus[i];
        new_emu.is_usefull = false;
        new_emu.network = network_num;
        all_known_emus.push_back(new_emu);
    }
}

int FindNetworkSize(const std::vector<EmuReaction>& reactions) {
    int current_size = 0;
    for (const bool state : reactions[0].right.emu.atom_states) {
        current_size += static_cast<int>(state);
    }
    return current_size;
}

Matrix GenerateDiffFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                              int rows, int cols, int id, int position, const std::vector<int>& id_to_pos,
                              const Matrix& nullspace) {
    Matrix matrix = Matrix::Zero(rows, cols);
    for (const FluxCombination& combination : symbolic_matrix) {
        double value = 0.0;
        for (const FluxAndCoefficient& flux : combination.fluxes) {
            if (id_to_pos[flux.id] == -1) {
                if (flux.id == id) {
                    value += flux.coefficient;
                }
            } else {
                value += -nullspace(id_to_pos[flux.id], position) * flux.coefficient;
            }

        }
        matrix(combination.i, combination.j) = value;
    }

    return matrix;
}
}
}
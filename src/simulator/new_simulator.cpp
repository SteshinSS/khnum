#include "simulator/new_simulator.h"

#include "utilities/debug_utills/debug_prints.h"

#include <iostream>
#include <tuple>
#include <set>


namespace khnum {

NewSimulator::NewSimulator(const std::vector<EmuNetwork> &networks, const std::vector<EmuAndMid> &input_mids,
                           const std::vector<Emu> &measured_isotopes) :
                    networks_{networks},
                    input_mids_{input_mids},
                    measured_isotopes_{measured_isotopes} {
    usefull_emus_.resize(networks_.size());
    unknown_size_.resize(networks_.size());
    known_size_.resize(networks_.size());
    network_size_.resize(networks_.size());
    mids_Yi_.resize(networks_.size());
    symbolic_Ai_.resize(networks_.size());
    symbolic_Bi_.resize(networks_.size());
    final_emus_.resize(networks_.size());

    std::vector<NetworkEmu> all_known_emus;
    for (size_t i = 0; i < input_mids_.size(); ++i) {
        NetworkEmu input_emu;
        input_emu.emu = input_mids_[i].emu;
        input_emu.network = -1;
        input_emu.order_in_usefull_emus = i;
        input_emu.order_in_X = i;
        input_emu.is_usefull = true;
        all_known_emus.push_back(input_emu);
    }

    for (network_ = 0; network_ < networks.size(); ++network_) {
        std::vector<Emu> unknown_emus;
        std::vector<Emu> known_emus;
        std::vector<Convolution> convolutions;
        FillEmuLists(unknown_emus, known_emus, convolutions, all_known_emus);
        CreateSymbolicMatrices(unknown_emus, known_emus, convolutions);
        FillFinalEmu(unknown_emus);
        InsertIntoAllKnownEmus(unknown_emus, all_known_emus);
        unknown_size_[network_] = unknown_emus.size();
        known_size_[network_] = known_emus.size();
        convolutions_.emplace_back(convolutions);
        network_size_[network_] = FindNetworkSize();
    }
}


void NewSimulator::FillEmuLists(std::vector<Emu> &unknown_emus, std::vector<Emu> &known_emus,
                                std::vector<Convolution> &convolutions, std::vector<NetworkEmu> &all_known_emus) {
    std::set<Emu> seen_emus;
    for (const EmuReaction &reaction : networks_[network_]) {
        if (reaction.left.size() == 1) {
            if (seen_emus.find(reaction.left[0].emu) == seen_emus.end()) {
                CheckAndInsertEmu(reaction.left[0].emu, all_known_emus, known_emus, unknown_emus);
                seen_emus.insert(reaction.left[0].emu);
            }
        } else {
            ConvolveReaction(reaction, all_known_emus, seen_emus, convolutions);
        }

        if (seen_emus.find(reaction.right.emu) == seen_emus.end()) {
            unknown_emus.push_back(reaction.right.emu);
            seen_emus.insert(reaction.right.emu);
        }
    }
}


void NewSimulator::CheckAndInsertEmu(const Emu &emu, std::vector<NetworkEmu> &all_known_emus,
                       std::vector<Emu> &known_emus, std::vector<Emu> &unknown_emus) {
    auto it = std::find_if(all_known_emus.begin(),
                           all_known_emus.end(),
                           [&emu](const NetworkEmu& known_emu) {
                               return known_emu.emu == emu;
                           });

    if (it != all_known_emus.end()) {
        if (!it->is_usefull) {
            it->is_usefull = true;
            usefull_emus_[it->network].push_back(it->order_in_X);
            it->order_in_usefull_emus = usefull_emus_[it->network].size() - 1;
        }
        known_emus.push_back(it->emu);
        mids_Yi_[network_].push_back({it->network, it->order_in_usefull_emus});
    } else {
        unknown_emus.push_back(emu);
    }
}


void NewSimulator::ConvolveReaction(const EmuReaction& reaction, std::vector<NetworkEmu> &all_known_emus,
                                           std::set<Emu>& seen_emus, std::vector<Convolution>& convolutions) {
    int reaction_id = reaction.id;

    auto it = std::find_if(convolutions.begin(), convolutions.end(),
                            [reaction_id](const Convolution& convolution) {
                                return convolution.flux_id == reaction_id;
                            });
    if (it != convolutions.end()) {
        return;
    }

    Convolution convolution;
    convolution.flux_id = reaction.id;
    for (const EmuSubstrate& emu : reaction.left) {
        if (seen_emus.find(emu.emu) != seen_emus.end()) {
            continue;
        }
        seen_emus.insert(emu.emu);

        auto it = std::find_if(all_known_emus.begin(),
                               all_known_emus.end(),
                               [&emu](const NetworkEmu& known_emu) {
                                   return known_emu.emu == emu.emu;
                               });

        if (!it->is_usefull) {
            it->is_usefull = true;
            usefull_emus_[it->network].push_back(it->order_in_X);
            it->order_in_usefull_emus = usefull_emus_[it->network].size() - 1;
        }

        PositionOfKnownEmu new_known_emu;
        new_known_emu.network = it->network;
        new_known_emu.position = it->order_in_usefull_emus;
        convolution.elements.emplace_back(new_known_emu);
    }

    convolutions.push_back(convolution);
}


void NewSimulator::CreateSymbolicMatrices(const std::vector<Emu>& unknown_emus,
                            const std::vector<Emu>& known_emus,
                            const std::vector<Convolution>& convolutions) {
    const int unknown_size = unknown_emus.size();
    const int known_size = known_emus.size() + convolutions.size();

    // A[unknown_size][unknown_size]
    std::vector<std::vector<FluxCombination>> A(unknown_size, std::vector<FluxCombination>(unknown_size));

    // B[unknown_size][known_size]
    std::vector<std::vector<FluxCombination>> B(unknown_size, std::vector<FluxCombination>(known_size));

    for (const EmuReaction &reaction : networks_[network_]) {
        int position_of_product = FindUnknownEmuPosition(reaction.right.emu, unknown_emus);

        {
            FluxAndCoefficient product;
            product.coefficient = -reaction.right.coefficient;
            product.id = reaction.id;
            A[position_of_product][position_of_product].fluxes.emplace_back(product);
        }

        if (reaction.left.size() > 1) {
            int position_of_convolution = FindConvolutionPosition(reaction.id, convolutions);
            FluxAndCoefficient convolution_element;
            convolution_element.coefficient = -reaction.right.coefficient;
            convolution_element.id = reaction.id;
            B[position_of_product][known_emus.size() + position_of_convolution].fluxes.emplace_back(convolution_element);
        } else {
            EmuSubstrate substrate = reaction.left[0];

            // returns -1 if emu is unknown
            int position_of_substrate = FindKnownEmuPosition(substrate.emu, known_emus);
            if (position_of_substrate == -1) {
                FluxAndCoefficient substrate_element;
                substrate_element.coefficient = reaction.right.coefficient;
                substrate_element.id = reaction.id;
                position_of_substrate = FindUnknownEmuPosition(substrate.emu, unknown_emus);
                A[position_of_product][position_of_substrate].fluxes.emplace_back(substrate_element);
            } else {
                FluxAndCoefficient substrate_element;
                substrate_element.coefficient = -substrate.coefficient;
                substrate_element.id = reaction.id;
                B[position_of_product][position_of_substrate].fluxes.emplace_back(substrate_element);
            }
        }
    }

    ConvertToSparseMatrix(A, symbolic_Ai_[network_]);
    ConvertToSparseMatrix(B, symbolic_Bi_[network_]);
}



int NewSimulator::FindUnknownEmuPosition(const Emu &emu,
                                         const std::vector<Emu>& unknown_emus) {
    auto position = find(unknown_emus.begin(),
                         unknown_emus.end(),
                         emu);

    return position - unknown_emus.begin();
}


int NewSimulator::FindKnownEmuPosition(const Emu &emu,
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


int NewSimulator::FindConvolutionPosition(const int reaction_id,
                                          const std::vector<Convolution>& convolutions) {
    auto position = find_if(convolutions.begin(),
                            convolutions.end(),
                            [reaction_id](const Convolution& convolution) {
                                return reaction_id == convolution.flux_id;
                            });

    return position - convolutions.begin();
}


void NewSimulator::ConvertToSparseMatrix(const std::vector<std::vector<FluxCombination>>& dense_matrix,
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


void NewSimulator::FillFinalEmu(const std::vector<Emu>& unknown_emus) {
    for (size_t i = 0; i < unknown_emus.size(); ++i) {
        Emu emu = unknown_emus[i];
        auto it = std::find_if(measured_isotopes_.begin(),
                               measured_isotopes_.end(),
                               [&emu](const Emu& measured_emu) {
                                   return measured_emu == emu;
                               });

        if (it != measured_isotopes_.end()) {
            FinalEmu final_emu;
            final_emu.emu = emu;
            final_emu.network = network_;
            final_emu.order_in_X = i;
            final_emu.position_in_result = it - measured_isotopes_.begin();
            final_emus_[network_].push_back(final_emu);
        }
    }
}


void NewSimulator::InsertIntoAllKnownEmus(std::vector<Emu>& unknown_emus,
                                          std::vector<NetworkEmu> &all_known_emus) {
    for (size_t i = 0; i < unknown_emus.size(); ++i) {
        NetworkEmu new_emu;
        new_emu.order_in_X = i;
        new_emu.emu = unknown_emus[i];
        new_emu.is_usefull = false;
        new_emu.network = network_;
        all_known_emus.push_back(new_emu);
    }
}


int NewSimulator::FindNetworkSize() {
    int current_size = 0;
    for (const bool state : networks_[network_][0].right.emu.atom_states) {
        current_size += static_cast<int>(state);
    }
    return current_size;
}


std::vector<EmuAndMid> NewSimulator::CalculateMids(const std::vector<Flux>& fluxes) {
    std::vector<EmuAndMid> result(measured_isotopes_.size());
    std::vector<std::vector<Mid>> known_mids(networks_.size());
    for (network_ = 0; network_ < networks_.size(); ++network_) {
        Matrix A = GenerateFluxMatrix(symbolic_Ai_[network_], fluxes, unknown_size_[network_]);
        Matrix B = GenerateFluxMatrix(symbolic_Bi_[network_], fluxes, known_size_[network_] + convolutions_[network_].size());
        Matrix Y = GenerateYMatrix(known_mids);
        Matrix BY = B * Y;
        // Matrix X = A.lu().solve(BY); // need to add -fopenmp to flags
        Matrix X = A.householderQr().solve(BY);
        // Matrix X = A.colPivHouseholderQr().solve(BY);
        SaveNewEmus(X, known_mids, result);
    }
    return result;
}


Matrix NewSimulator::GenerateFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                                        const std::vector<Flux>& fluxes,
                                        const int cols) {
    Matrix matrix =  Matrix::Zero(unknown_size_[network_], cols);
    for (const FluxCombination& combination : symbolic_matrix) {
        double value = 0.0;
        for (const FluxAndCoefficient& flux : combination.fluxes) {
            value += flux.coefficient * fluxes[flux.id];
        }
        matrix(combination.i, combination.j) = value;
    }

    return matrix;
}


Matrix NewSimulator::GenerateYMatrix(const std::vector<std::vector<Mid>>& known_mids) {
    Matrix Y(known_size_[network_] + convolutions_[network_].size(), network_size_[network_] + 1);
    for (size_t i = 0; i < mids_Yi_[network_].size(); ++i) {
        PositionOfKnownEmu known_emu = mids_Yi_[network_][i];
        Mid mid;
        if (known_emu.network == -1) {
            mid = input_mids_[known_emu.position].mid;
        } else {
            mid = known_mids[known_emu.network][known_emu.position];
        }
        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y(i, mass_shift) = mid[mass_shift];
        }
    }

    for (size_t i = 0; i < convolutions_[network_].size(); ++i) {
        Mid mid(1, 1.0); // MID = [1.0]
        for (const PositionOfKnownEmu& emu : convolutions_[network_][i].elements) {
            if (emu.network == -1) {
                mid = mid * input_mids_[emu.position].mid;
            } else {
                mid = mid * known_mids[emu.network][emu.position];
            }
        }
        for (size_t mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y(i + known_size_[network_], mass_shift) = mid[mass_shift];
        }
    }

    return Y;
}


void NewSimulator::SaveNewEmus(const Matrix& X,
                               std::vector<std::vector<Mid>>& known_mids,
                               std::vector<EmuAndMid>& result) {

    for (int position : usefull_emus_[network_]) {
        Mid new_mid;
        for (int mass_shift = 0; mass_shift < network_size_[network_] + 1; ++mass_shift) {
            new_mid.push_back(X(position, mass_shift));
        }
        known_mids[network_].push_back(new_mid);
    }

    for (FinalEmu& final_emu : final_emus_[network_]) {
        Mid result_mid;
        for (int mass_shift = 0; mass_shift < network_size_[network_] + 1; ++mass_shift) {
            result_mid.push_back(X(final_emu.order_in_X, mass_shift));
        }
        EmuAndMid result_emu;
        result_emu.emu = final_emu.emu;
        result_emu.mid = result_mid;
        result[final_emu.position_in_result] = result_emu;
    }
}


} // namespace khnum

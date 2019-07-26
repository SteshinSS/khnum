#include "new_simulator.h"

#include <iostream>


namespace khnum {

NewSimulator::NewSimulator(const std::vector<EmuNetwork> &networks, const std::vector<EmuAndMid> &input_mids,
                           const std::vector<Emu> &measured_isotopes) :
                    networks_{networks},
                    input_mids_{input_mids},
                    measured_isotopes_{measured_isotopes} {
    usefull_emus_.resize(networks_.size());
    unknown_size_.resize(networks_.size());
    known_size_.resize(networks_.size());

    std::vector<NetworkEmu> all_known_emus;
    for (network_ = 0; network_ < networks.size(); ++network_) {
        std::vector<Emu> unknown_emus;
        std::vector<Emu> known_emus;
        std::vector<Convolution> convolutions;
        FillEmuLists(unknown_emus, known_emus, convolutions, all_known_emus);
        CreateSymbolicMatrices(unknown_emus, known_emus, convolutions);
        unknown_size_[network_] = unknown_emus.size();
        known_size_[network_] = known_emus.size();
        convolutions_.emplace_back(convolutions);
    }
}


void NewSimulator::FillEmuLists(std::vector<Emu> &unknown_emus, std::vector<Emu> &known_emus,
                                std::vector<Convolution> &convolutions, std::vector<NetworkEmu> &all_known_emus) {
    // Fills known_emus and unknown_emus
    for (const EmuReaction &reaction : networks_[network_]) {
        if (reaction.left.size() == 1) {
            CheckAndInsertEmu(reaction.left[0].emu, all_known_emus, known_emus, unknown_emus);
        } else {
            Convolution new_convolution = ConvolveReaction(reaction, all_known_emus);
            convolutions.push_back(new_convolution);
        }

        unknown_emus.push_back(reaction.right.emu);
    }

    // delete repeated emus
    std::sort(known_emus.begin(), known_emus.end());
    known_emus.erase(std::unique(known_emus.begin(), known_emus.end()), known_emus.end());

    std::sort(unknown_emus.begin(), unknown_emus.end());
    unknown_emus.erase(std::unique(unknown_emus.begin(), unknown_emus.end()), unknown_emus.end());

    std::sort(known_emus_[network_].begin(), known_emus_[network_].end(),
                [](const PositionOfKnownEmu& lhs, const PositionOfKnownEmu& rhs) {
                    return std::tie(lhs.network, lhs.position) < std::tie(rhs.network, rhs.position);
                });
    auto known_last = std::unique(known_emus_.begin(), known_emus_.end(),
                [](const PositionOfKnownEmu& lhs, const PositionOfKnownEmu& rhs) {
                    return std::tie(lhs.network, lhs.position) == std::tie(rhs.network, rhs.position);
                });

    known_emus_.erase(known_last, known_emus_.end());


    std::sort(convolutions.begin(), convolutions.end(),
                [](const Convolution& lhs, const Convolution& rhs) {
                    return lhs.flux_id < rhs.flux_id;
                });

    auto last = std::unique(convolutions.begin(), convolutions.end(),
                            [](const Convolution& lhs, const Convolution& rhs) {
                                return lhs.flux_id == rhs.flux_id;
                            });

    convolutions.erase(last, convolutions.end());
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
            it->order_in_usefull_emus = usefull_emus_.size() - 1;
        }
        known_emus.push_back(it->emu);
        known_emus_[network_].push_back({it->network, it->order_in_usefull_emus});
    } else {
        unknown_emus.push_back(it->emu);
    }
}


Convolution NewSimulator::ConvolveReaction(const EmuReaction& reaction, std::vector<NetworkEmu> &all_known_emus) {
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
            usefull_emus_[it->network].push_back(it->order_in_X);
            it->order_in_usefull_emus = usefull_emus_.size() - 1;
        }

        PositionOfKnownEmu new_known_emu;
        new_known_emu.network = it->network;
        new_known_emu.position = it->order_in_usefull_emus;
        convolution.elements.emplace_back(new_known_emu);
    }

    return convolution;
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
        FluxAndCoefficient product;
        product.coefficient = -reaction.right.coefficient;
        product.id = reaction.id;
        A[position_of_product][position_of_product].fluxes.emplace_back(product);

        if (reaction.left.size() > 1) {
            int position_of_convolution = FindConvolutionPosition(reaction.id, convolutions);
            FluxAndCoefficient convolution_element;
            convolution_element.coefficient = 1.0;
            convolution_element.id = reaction.id;
            B[position_of_product][known_emus.size() + position_of_convolution].fluxes.emplace_back(convolution_element);
        } else {
            EmuSubstrate substrate = reaction.left[0];

            // return -1 if emu is unknown
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
                            [&emu](const EmuAndMid &known_mid) {
                                return known_mid.emu == emu;
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
    for (int i = 0; i < dense_matrix.size(); ++i) {
        for (int j = 0; j < dense_matrix[0].size(); ++j) {
            if (!dense_matrix[i][j].fluxes.empty()) {
                FluxCombination matrix_element = dense_matrix[i][j];
                matrix_element.i = i;
                matrix_element.j = j;
                sparse_matrix.emplace_back(matrix_element);
            }
        }
    }
}

std::vector<EmuAndMid> NewSimulator::CalculateMids(const std::vector<Flux>& fluxes) {
    std::vector<std::vector<Mid>> known_mids;
    for (network_ = 0; network_ < networks_.size(); ++network_) {
        Matrix A = GenerateFluxMatrix(symbolic_Ai_[network_], fluxes);
        Matrix B = GenerateFluxMatrix(symbolic_Bi_[network_], fluxes);
        Matrix Y = GenerateYMatrix(known_mids);
        Matrix BY = B * Y;
        Matrix X = A.householderQr().solve(BY);
        SaveNewEmus(X, known_mids);
    }
}


Matrix NewSimulator::GenerateFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                                        const std::vector<Flux>& fluxes) {
    Matrix matrix =  Matrix::Zero(unknown_emus.size(), unknown_emus.size());
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
    Matrix Y(a, a);
    for (int i = 0; i < symbolic_Yi_[network_].size(); ++i) {
        const YElement& new_component = symbolic_Yi_[network_][i];
        const Mid& mid = known_mids[new_component.network][new_component.position];
        for (int mass_shift = 0; mass_shift < mid.size(); ++mass_shift) {
            Y(i, mass_shift) = mid[mass_shift];
        }
    }

    return Y;
}


void NewSimulator::SaveNewEmus(const Matrix& X,
                 std::vector<std::vector<Mid>>& known_mids) {
    for (int i : usefull_emus_[network_]) {
        Mid usefull_mid(X.cols());
        for (int mass_shift = 0; mass_shift < X.cols(); ++mass_shift) {
            usefull_mid.push_back(X(i, mass_shift));
        }
        known_mids[network_].push_back(usefull_mid);
    }

    for (const YConvolution& convolution : convolutions_[network_]) {

    }
}


} // namespace khnum

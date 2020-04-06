#pragma once

#include <vector>
#include <unordered_map>

#include "utilities/matrix.h"
#include "utilities/emu.h"
#include "simulator/flux_combination.h"
#include "simulator/simulation_data.h"
#include "utilities/measurement.h"


namespace khnum {
namespace generator_utilites {
void FillEmuLists(const std::vector<EmuReaction>& reactions,
                  std::vector<NetworkEmu>& all_known_emus,
                  GeneratorNetworkData& network_data,
                  std::vector<std::vector<int>>& usefull_emus);

void CheckAndInsertEmu(const Emu &emu,
                       std::vector<NetworkEmu>& all_known_emus,
                       GeneratorNetworkData& network_data,
                       std::vector<std::vector<int>>& usefull_emus);


Convolution ConvolveReaction(const EmuReaction& reaction,
                             std::vector<NetworkEmu>& all_known_emus,
                             std::vector<std::vector<int>>& usefull_emus);

void CreateSymbolicMatrices(const std::vector<EmuReaction>& reactions,
                            GeneratorNetworkData& network_data);

int FindUnknownEmuPosition(const Emu &emu,
                           const std::vector<Emu>& unknown_emus);

int FindKnownEmuPosition(const Emu &emu,
                         const std::vector<Emu>& known_emus);

void ConvertToSparseMatrix(const std::vector<std::vector<FluxCombination>>& dense_matrix,
                           std::vector<FluxCombination>& sparse_matrix);

void FillFinalEmu(const std::vector<Measurement>& measured_isotopes,
                  GeneratorNetworkData& network_data);

void InsertIntoAllKnownEmus(const std::vector<Emu>& unknown_emus,
                            int network_num,
                            std::vector<NetworkEmu>& all_known_emus);

int FindNetworkSize(const std::vector<EmuReaction>& reactions);

Matrix GenerateDiffFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                              int rows, int cols, int id, int position, const std::vector<int>& id_to_pos,
                              const Matrix& nullspace);
}
}
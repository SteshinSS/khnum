#pragma once

#include <vector>
#include <set>

#include "simulator/flux_combination.h"
#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"
#include "utilities/problem.h"


/*
namespace khnum {

class OldSimulator {
public:
    OldSimulator(const SimulatorParameters& parameters);

    void CalculateMids(const std::vector<Flux> &fluxes);

private:
    // Preprocessing functions
    void FillEmuLists(std::vector<Emu> &unknown_emus, std::vector<Emu> &known_emus,
                      std::vector<Convolution> &convolutions);

    void CheckAndInsertEmu(const Emu &emu, std::vector<Emu> &known_emus, std::vector<Emu> &unknown_emus);

    void InsertIntoAllKnownEmus(std::vector<Emu>& unknown_emus);

    void FillFinalEmu(const std::vector<Emu>& unknown_emus);

    Convolution ConvolveReaction(const EmuReaction& reaction);

    void CreateSymbolicMatrices(const std::vector<Emu>& unknown_emus,
                                const std::vector<Emu>& known_emus,
                                const std::vector<Convolution>& convolutions);

    int FindUnknownEmuPosition(const Emu &emu,
                               const std::vector<Emu>& unknown_emus);

    int FindKnownEmuPosition(const Emu &emu,
                             const std::vector<Emu>& known_emus);

    void ConvertToSparseMatrix(const std::vector<std::vector<FluxCombination>>& dense_matrix,
                               std::vector<FluxCombination>& sparse_matrix);

    int FindNetworkSize();

    // Calculations functions
    Matrix GenerateFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                              const std::vector<Flux>& fluxes,
                              const int cols);

    Matrix GenerateYMatrix(const std::vector<std::vector<Mid>> &known_mids);

    void SaveNewEmus(const Matrix &X, std::vector<std::vector<Mid>> &known_mids, std::vector<EmuAndMid>& result);


    Matrix GenerateDiffFluxMatrix(const std::vector<FluxCombination>& symbolic_matrix,
                                                int cols,
                                                int id, int position);

    Matrix GenerateDiffYMatrix(const std::vector<std::vector<Mid>>& known_d_mids,
                                             std::vector<std::vector<Mid>>& known_mids);

private:
    size_t network_;
    Matrix nullspace_;
    const std::vector<int> id_to_pos_;
    const std::vector<int>& free_fluxes_id_;
    bool use_analytic_jacobian_;


private:
    const std::vector<EmuNetwork> networks_;
    const std::vector<EmuAndMid> input_mids_;
    const std::vector<Emu> measured_isotopes_;
    std::vector<NetworkEmu> all_known_emus_;

    // number of known/unknown emus in the i'th network
    std::vector<int> unknown_size_; //
    std::vector<int> known_size_; //
    std::vector<int> network_size_;

    std::vector<int> depended_flux_positions_;

    std::vector<std::vector<FluxCombination>> symbolic_Ai_; //
    std::vector<std::vector<FluxCombination>> symbolic_Bi_; //

    // This forms Y matrix
    std::vector<std::vector<PositionOfSavedEmu>> mids_Yi_; //
    std::vector<std::vector<Convolution>> convolutions_;  //

    // vector usefull_emu_[i] contains positions of Xi emus, which are using later
    std::vector<std::vector<int>> usefull_emus_; //

    std::vector<std::vector<FinalEmu>> final_emus_;

};
} // namespace khnum

*/
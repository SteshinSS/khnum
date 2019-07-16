#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"

class Simulator {
public:
    Simulator(const std::vector<Flux> &fluxes, const std::vector<EmuNetwork> &networks,
              const std::vector<EmuAndMid> &all_known_emus, const std::vector<Emu> &measured_isotopes);

    std::vector<EmuAndMid> CalculateMids();

private:

    static int FindNetworkSize(const EmuNetwork &network);
    std::vector<EmuAndMid> SelectMeasuredMID();

    void SolveOneNetwork(const EmuNetwork &network);


    void FillEmuLists(std::vector<Emu> &unknown_emus,
                      std::vector<EmuAndMid> &known_emus,
                      const EmuNetwork &network);

    void CheckIsEmuKnown(const Emu &emu, const std::vector<EmuAndMid> &where_find_emu_list,
                         std::vector<EmuAndMid> &known_emus, std::vector<Emu> &unknown_emus);


    static const Mid *FindMid(const Emu &emu,
                       const std::vector<EmuAndMid> &mids);


    EmuAndMid ConvolveEmu(const EmuReactionSide &convolve_reaction);

    static Matrix FormYMatrix(const std::vector<EmuAndMid> &known_emus,
                       const int current_size);

    void FillABMatrices(Matrix &A, Matrix &B,
                        const EmuNetwork &network,
                        const std::vector<EmuAndMid> &known_emus,
                        const std::vector<Emu> &unknown_emus);

    static int FindUnknownEmuPosition(const Emu &emu,
                               const std::vector<Emu> unknown_emus);


    static int FindKnownEmuPosition(const Emu &emu,
                             const std::vector<EmuAndMid> known_emus);

    void AppendNewMids(const Matrix &X,
                       const std::vector<Emu> &unknown_emus,
                       const int current_size);


private:
    const std::vector<Flux> fluxes_;
    const std::vector<EmuNetwork> networks_;
    std::vector<EmuAndMid> all_known_emus_;
    const std::vector<Emu> measured_isotopes_;
};

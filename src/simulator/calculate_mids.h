#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"


std::vector<EmuAndMid> CalculateMids(const std::vector<Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EmuAndMid> all_known_emus,
                                     const std::vector<Emu> &measured_isotopes);


int FindNetworkSize(const EMUNetwork &network);


std::vector<EmuAndMid> SelectMeasuredMID(const std::vector<EmuAndMid> &all_known_emus,
                                         const std::vector<Emu> &measured_isotopes);

void SolveOneNetwork(const std::vector<Flux> &fluxes,
                     const EMUNetwork &network,
                     std::vector<EmuAndMid> &all_known_emus);


void FillEMULists(std::vector<Emu> &unknown_emus,
                  std::vector<EmuAndMid> &known_emus,
                  const EMUNetwork &network,
                  const std::vector<EmuAndMid> &all_known_emus);


void CheckIsEmuKnown (const Emu& emu, const std::vector<EmuAndMid>& where_find_emu_list,
                      std::vector<EmuAndMid>& known_emus, std::vector<Emu>& unknown_emus);

const Mid *FindMid(const Emu &emu,
                   const std::vector<EmuAndMid> &known_mids);

EmuAndMid ConvolveEmu(const EmuReactionSide &convolve_reaction,
                      const std::vector<EmuAndMid> &known_mids);

Matrix FormYMatrix(const std::vector<EmuAndMid> &known_emus,
                   const int current_size);

void FillABMatrices(Matrix &A, Matrix &B,
                    const EMUNetwork &network,
                    const std::vector<EmuAndMid> &known_emus,
                    const std::vector<Emu> &unknown_emus,
                    const std::vector<Flux> &fluxes,
                    const std::vector<EmuAndMid> &known_mids);

int FindUnknownEmuPosition(const Emu &emu,
                           const std::vector<Emu> unknown_emus);


int FindKnownEmuPosition(const Emu &emu,
                         const std::vector<EmuAndMid> known_emus);

void AppendNewMids(const Matrix &X,
                   const std::vector<Emu> &unknown_emus,
                   std::vector<EmuAndMid> &known_mids,
                   const int current_size);


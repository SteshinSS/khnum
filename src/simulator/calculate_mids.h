#pragma once

#include <vector>

#include "utilities/emu.h"
#include "utilities/emu_and_mid.h"
#include "utilities/reaction.h"
#include "utilities/matrix.h"


std::vector<EmuAndMid> CalculateMids(const std::vector<Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EmuAndMid> known_mids,
                                     const std::vector<Emu> &measured_isotopes);


int FindNetworkSize(const EMUNetwork &network);


std::vector<EmuAndMid> SelectMeasuredMID(const std::vector<EmuAndMid> &known_mids,
                                         const std::vector<Emu> &measured_isotopes);

void SolveOneNetwork(const std::vector<Flux> &fluxes,
                     const EMUNetwork &network,
                     std::vector<EmuAndMid> &known_mids);


void FillEMULists(std::vector<Emu> &unknown_emus,
                  std::vector<EmuAndMid> &known_emus,
                  const EMUNetwork &network,
                  const std::vector<EmuAndMid> &known_mids);


void FormYMatrix(Matrix &Y,
                 const std::vector<EmuAndMid> &known_emus,
                 const int current_size);


void FormABMatrices(Matrix &A, Matrix &B,
                    const EMUNetwork &network,
                    const std::vector<EmuAndMid> &known_emus,
                    const std::vector<Emu> &unknown_emus,
                    const std::vector<Flux> &fluxes,
                    const std::vector<EmuAndMid> &known_mids);


void AppendNewMIDs(const Matrix &X,
                   const std::vector<Emu> &unknown_emus,
                   std::vector<EmuAndMid> &known_mids,
                   const int current_size);

bool IsEMUKnown(const Emu &emu,
                const std::vector<EmuAndMid> known_emus);


int FindUnknownEMUsPosition(const Emu &emu,
                            const std::vector<Emu> unknown_emus);


int FindKnownEMUsPosition(const Emu &emu,
                          const std::vector<EmuAndMid> known_emus);

const Mid *GetMID(const Emu &emu,
                  const std::vector<EmuAndMid> &known_mids);

EmuAndMid ConvolveEMU(const EmuReactionSide &convolve_reaction,
                      const std::vector<EmuAndMid> &known_mids);

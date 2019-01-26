#ifndef CFLEX_CALCULATE_MIDS_H
#define CFLEX_CALCULATE_MIDS_H


#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"
#include "../math/math_utilites.h"

#include <vector>
#include <map>

std::vector<EMUandMID> CalculateMids(const std::map<std::string, Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     std::vector<EMUandMID> known_mids,
                                     const std::vector<EMU> &measured_isotopes);


int FindNetworkSize(const EMUNetwork &network);


std::vector<EMUandMID> SelectMeasuredMID(const std::vector<EMUandMID> &known_mids,
                                         const std::vector<EMU> &measured_isotopes);

void SolveOneNetwork(const std::map<std::string, Flux> &fluxes,
                     const EMUNetwork &network,
                     std::vector<EMUandMID> &known_mids);


void FillEMULists(std::vector<EMU> &unknown_emus,
                  std::vector<EMUandMID> &known_emus,
                  const EMUNetwork &network,
                  const std::vector<EMUandMID> &known_mids);


void FormYMatrix(Matrix &Y,
                 const std::vector<EMUandMID> &known_emus,
                 const int current_size);


void FormABMatrices(Matrix &A, Matrix &B,
                    const EMUNetwork &network,
                    const std::vector<EMUandMID> &known_emus,
                    const std::vector<EMU> &unknown_emus,
                    const std::map<std::string, Flux> &fluxes,
                    const std::vector<EMUandMID> &known_mids);


void AppendNewMIDs(const Matrix &X,
                   const std::vector<EMU> &unknown_emus,
                   std::vector<EMUandMID> &known_mids,
                   const int current_size);

bool IsEMUKnown(const EMU &emu,
                const std::vector<EMUandMID> known_emus);


int FindUnknownEMUsPosition(const EMU &emu,
                            const std::vector<EMU> unknown_emus);


int FindKnownEMUsPosition(const EMU &emu,
                          const std::vector<EMUandMID> known_emus);

const MID *GetMID(const EMU &emu,
                  const std::vector<EMUandMID> &known_mids);

EMUandMID ConvolveEMU(const EMUReactionSide &convolve_reaction,
                      const std::vector<EMUandMID> &known_mids);

#endif //CFLEX_CALCULATE_MIDS_H

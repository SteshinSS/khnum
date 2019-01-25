#ifndef CFLEX_CALCULATE_MIDS_H
#define CFLEX_CALCULATE_MIDS_H


#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include "../utilities/reaction_struct.h"

#include <vector>

std::vector<EMUandMID> CalculateMids(const std::vector<Flux> &fluxes,
                                     const std::vector<EMUNetwork> &networks,
                                     const std::vector<EMUandMID> &input_substrates_mids,
                                     const std::vector<EMU> &measured_isotopes);


std::vector<EMUandMID> SolveOneNetwork(const std::vector<Flux> &fluxes,
                                       const EMUNetwork &network,
                                       const std::vector<EMUandMID> &known_mids,
                                       int current_size);

const MID * GetMID(const EMU &emu,
                const std::vector<EMUandMID> &known_mids);

EMUandMID ConvolveEMU(const EMUReactionSide &convolve_reaction,
                      const std::vector<EMUandMID> &known_mids);
#endif //CFLEX_CALCULATE_MIDS_H

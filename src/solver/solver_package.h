#ifndef CFLEX_SOLVER_PACKAGE_H
#define CFLEX_SOLVER_PACKAGE_H


struct SolverPackage {
    Matrix* const nullspace;
    std::vector<Reaction>* const reactions;
    std::vector<EMUNetwork>* const networks;
    std::vector<EmuAndMid>* const input_mids;
    std::vector<Emu>* const measured_isotopes;
    std::vector<Measurement>* const measurements;
};


#endif //CFLEX_SOLVER_PACKAGE_H

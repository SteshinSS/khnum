#ifndef CFLEX_MODELLER_H
#define CFLEX_MODELLER_H


#include "utilities/problem.h"

#include "parser_results.h"

class Modeller {
public:
    Modeller(const ParserResults);

    void CreateEmuNetworks();
    void CalculateInputSubstrateMids();
    void CreateNullspaceMatrix();
    void CalculateFluxBounds();
    
    Problem GetProblem();

private:
    std::vector<Reaction> reactions_;
    std::vector<Emu> measured_isotopes_;
    std::vector<Measurement> measurements_;
    std::vector<InputSubstrate> input_substrate_;
    std::vector<std::string> excluded_metabolites_;

    std::vector<Emu> input_emu_list_;
    std::vector<EMUReaction> all_emu_reactions_;

    Matrix nullspace_;
    std::vector<EmuAndMid> input_substrate_mids_;
    std::vector<EMUNetwork> emu_networks_;
};

#endif //CFLEX_MODELLER_H

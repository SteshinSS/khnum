#ifndef CFLEX_PARSER_RESULTS_H
#define CFLEX_PARSER_RESULTS_H

#include <vector>
#include <optional>

#include "reaction_struct.h"
#include "Emu.h"
#include "measurements_struct.h"
#include "input_substrate.h"

struct ParserResults {
    std::vector<Reaction> reactions;
    std::vector<Emu> measured_isotopes;
    std::vector<Measurement> measurements;
    std::vector<InputSubstrate> input_substrate;
    std::vector<std::string> excluded_metabolites;
};

#endif //CFLEX_PARSER_RESULTS_H

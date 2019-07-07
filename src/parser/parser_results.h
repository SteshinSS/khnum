#ifndef CFLEX_PARSER_RESULTS_H
#define CFLEX_PARSER_RESULTS_H

#include <vector>
#include <optional>

#include "reaction_struct.h"
#include "Emu.h"
#include "measurements_struct.h"
#include "input_substrate.h"

struct ParserResults {
    std::optional<std::vector<Reaction>> reactions;
    std::optional<std::vector<Emu>> measuredEmu;
    std::optional<std::vector<Measurement>> measurements;
    std::optional<std::vector<InputSubstrate>> input_substrate;
    std::optional<std::vector<std::string>> excluded_metabolites;
};

#endif //CFLEX_PARSER_RESULTS_H

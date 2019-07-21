#pragma once

#include <vector>
#include <optional>

#include "utilities/reaction.h"
#include "utilities/emu.h"
#include "utilities/measurement.h"
#include "utilities/input_substrate.h"


namespace khnum {
struct ParserResults {
    std::vector<Reaction> reactions;
    std::vector<Emu> measured_isotopes;
    std::vector<Measurement> measurements;
    std::vector<InputSubstrate> input_substrate;
    std::vector<std::string> excluded_metabolites;
};
} //namespace khnum

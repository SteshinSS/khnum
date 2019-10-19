#pragma once

#include <vector>
#include <string>
#include <utilities/reaction.h>

#include "utilities/measurement.h"
#include "utilities/input_substrate.h"
#include "utilities/emu.h"


namespace khnum {
namespace open_flux_parser {

struct Delimiters {
    char csv_delimiter;
    std::string reaction_side_delimiter;
    std::string substrate_delimiter;
};

std::vector<std::string> GetLines(const std::string &path);

std::vector<Emu> ParseMeasuredIsotopes(const std::vector<std::string> &raw_measured_isotopes);

Emu ParseOneMeasuredIsotope(const std::string& raw_measured_isotope, int line_number = -1 );

std::vector<Measurement> ParseMeasurements(const std::vector<std::string>& raw_measurements,
                                           const std::vector<Emu>& measured_isotopes,
                                           const Delimiters& delimiters);

std::string GetCell(std::stringstream &line, const Delimiters& delimiters);

std::vector<InputSubstrate> ParseInputSubstrates(const std::vector<std::string>& raw_input_substrates,
                                                 const Delimiters& delimiters);



std::vector<Reaction> ParseReactions(const std::string& reactions_path,
                                     const Delimiters& delimiters);

Reaction FillReaction(const std::string& raw_line, int id, const Delimiters& delimiters);

ChemicalEquation ParseChemicalEquation(std::stringstream &line, const Delimiters& delimiters);

ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation,
                                      const Delimiters& delimiters);

ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation, const Delimiters& delimiters);

void ParseAtomEquationSide(const std::string &atom_equation, const Delimiters& delimiters,
                           ChemicalEquationSide *equation_side);

ReactionType ParseReactionType(const std::string &type);

std::optional<Basis> ParseBasis(const std::string &basis);

Deviation ParseDeviation(const std::string &deviation);

} // open_flux_parser
} // khnum
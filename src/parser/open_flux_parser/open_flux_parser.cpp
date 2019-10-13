#include "parser/open_flux_parser/open_flux_parser.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>
#include <sstream>
#include <tuple>

#include "parser/open_flux_parser/open_flux_utills.h"
#include "modeller/sort_reactions.h"


namespace khnum {

ParserResults ParserOpenFlux::GetResults() {
    ParserResults results;
    results.reactions = reactions_;
    results.measured_isotopes = measured_isotopes_;
    results.measurements = measurements_;
    results.input_substrate = input_substrates_;
    results.excluded_metabolites = excluded_metabolites_;

    return results;
}


void ParserOpenFlux::ParseExcludedMetabolites() {
    const std::string excluded_metabolites_path = path_ + "/excluded_metabolites.txt";
    excluded_metabolites_ = open_flux_parser::GetLines(excluded_metabolites_path);
}


void ParserOpenFlux::ParseMeasuredIsotopes() {
    const std::string measured_isotopes_path = path_ + "/measured_isotopes.txt";

    std::vector<std::string> raw_measured_isotopes = open_flux_parser::GetLines(measured_isotopes_path);
    measured_isotopes_ = open_flux_parser::ParseMeasuredIsotopes(raw_measured_isotopes);
}


void ParserOpenFlux::ParseMeasurements() {
    const std::string measurements_path = path_ + "/measurements.csv";

    std::vector<std::string> raw_measurements = open_flux_parser::GetLines(measurements_path);
    measurements_ = open_flux_parser::ParseMeasurements(raw_measurements, measured_isotopes_, delimiters_);
}


void ParserOpenFlux::ParseSubstrateInput() {
    const std::string input_substrates_path = path_ + "/substrate_input.csv";
    const std::vector<std::string>& raw_substrates = open_flux_parser::GetLines(input_substrates_path);
    input_substrates_ = open_flux_parser::ParseInputSubstrates(raw_substrates, delimiters_);
}

void ParserOpenFlux::ParseReactions() {
    const std::string reactions_path = path_ + "/model.csv";
    reactions_ = open_flux_parser::ParseReactions(reactions_path, delimiters_);
    reactions_ = modelling_utills::SortReactionsByType(reactions_);
}


} //namespace khnum
#include "parser/open_flux_parser/open_flux_parser.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <utilities/debug_utills/debug_prints.h>

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


void ParserOpenFlux::Parse() {

    ParseMeasuredIsotopes();
    ParseMeasurements();
    ParseCorrectionMatrices();
    ParseSubstrateInput();
    ParseReactions();
    ParseExcludedMetabolites();
}


void ParserOpenFlux::ParseExcludedMetabolites() {
    const std::string excluded_metabolites_path = path_ + "/excluded_metabolites.txt";
    excluded_metabolites_ = open_flux_parser::GetLines(excluded_metabolites_path);

    std::unordered_set<std::string> left;
    std::unordered_set<std::string> right;

    for (const Reaction& reaction : reactions_) {
        for (const Substrate& substrate : reaction.chemical_equation.left) {
            if (std::find(excluded_metabolites_.begin(), excluded_metabolites_.end(), substrate.name) == excluded_metabolites_.end()) {
                left.insert(substrate.name);
            }

        }
        for (const Substrate& substrate : reaction.chemical_equation.right) {
            if (std::find(excluded_metabolites_.begin(), excluded_metabolites_.end(), substrate.name) == excluded_metabolites_.end()) {
                right.insert(substrate.name);
            }
        }
    }

    for (const std::string& substrate : left) {
        if (right.find(substrate) == right.end()) {
            excluded_metabolites_.push_back(substrate);
        }
    }
    for (const std::string& substrate : right) {
        if (left.find(substrate) == left.end()) {
            excluded_metabolites_.push_back(substrate);
        }
    }

    std::sort(excluded_metabolites_.begin(), excluded_metabolites_.end());
    excluded_metabolites_.erase(unique(excluded_metabolites_.begin(), excluded_metabolites_.end()), excluded_metabolites_.end());
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


void ParserOpenFlux::ParseCorrectionMatrices() {
    const std::string correction_matrices_path = path_ + "/correctionMatrices/";
    for (Measurement& measurement : measurements_) {
        std::string file = correction_matrices_path + open_flux_parser::GetMeasuredIsotopeName(measurement.emu) + ".txt";
        try {
            std::vector<std::string> raw_matrix = open_flux_parser::GetLines(file);
            measurement.correction_matrix = open_flux_parser::ParseCorrectionMatrix(raw_matrix, delimiters_);
        } catch (std::runtime_error& error) {
            continue;
        }
    }
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
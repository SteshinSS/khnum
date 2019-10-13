#include "parser/open_flux_parser.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>
#include <sstream>
#include <tuple>

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
    excluded_metabolites_ = ReadEachLine(excluded_metabolites_path);
}


void ParserOpenFlux::ParseMeasuredIsotopes() {
    const std::string measured_isotopes_path = path_ + "/measured_isotopes.txt";

    std::vector<std::string> raw_measured_isotopes = ReadEachLine(measured_isotopes_path);

    // line_number need for error messages
    int line_number = 1;
    for (const std::string &raw_measured_isotope : raw_measured_isotopes) {
        std::stringstream line(raw_measured_isotope);
        Emu new_emu;
        getline(line, new_emu.name, ':');

        std::string row_atom_states;
        getline(line, row_atom_states);

        new_emu.atom_states.resize(row_atom_states.size());
        for (size_t atom_position = 0; atom_position < row_atom_states.size(); ++atom_position) {
            if (row_atom_states[atom_position] == '1') {
                new_emu.atom_states[atom_position] = true;
            } else if (row_atom_states[atom_position] == '0') {
                new_emu.atom_states[atom_position] = false;
            } else {
                throw std::runtime_error(
                    "Line: " + std::to_string(line_number) + "There is strange atom states in measured isotope!");
            }
        }
        ++line_number;

        measured_isotopes_.push_back(new_emu);
    }
}


void ParserOpenFlux::ParseMeasurements() {
    const std::string measurements_path = path_ + "/measurements.csv";

    std::ifstream input(measurements_path);
    // skip table's head line
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    for (const Emu &isotope : measured_isotopes_) {
        Measurement new_measurement;
        new_measurement.emu = isotope;
        for (size_t mass_shift = 0; mass_shift < isotope.atom_states.size() + 1; ++mass_shift) {
            getline(input, raw_line);
            std::stringstream line(raw_line);

            double value = std::stod(GetCell(line));
            new_measurement.mid.push_back(value);

            double error = std::stod(GetCell(line));
            new_measurement.errors.push_back(error);
        }
        measurements_.push_back(new_measurement);
    }
}


void ParserOpenFlux::ParseSubstrateInput() {
    const std::string input_substrates_path = path_ + "/substrate_input.csv";

    std::ifstream input(input_substrates_path);

    // skip table head-line
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    while (getline(input, raw_line)) {
        std::stringstream line(raw_line);
        std::string input_substrate_name;
        getline(line, input_substrate_name, csv_delimiter_);

        auto input_substrate_iterator = std::find_if(input_substrates_.begin(),
                                                     input_substrates_.end(),
                                                     [&input_substrate_name](InputSubstrate &input_substrate) {
                                                         return input_substrate.name == input_substrate_name;
                                                     });

        if (input_substrate_iterator == input_substrates_.end()) {
            InputSubstrate new_input_substrate;
            new_input_substrate.name = input_substrate_name;
            input_substrates_.push_back(new_input_substrate);

            //move iterator back, so we can work with the new substrate
            input_substrate_iterator = input_substrates_.end();
            --input_substrate_iterator;
        }

        std::string raw_labeling_pattern;
        getline(line, raw_labeling_pattern, csv_delimiter_);
        std::stringstream labeling_pattern(raw_labeling_pattern);

        std::vector<Fraction> fractions;
        Fraction new_fraction;
        while (labeling_pattern >> new_fraction) {
            fractions.push_back(new_fraction);
        }

        std::string raw_ratio;
        getline(line, raw_ratio);
        double ratio = std::stod(raw_ratio);

        Mixture new_mixture;
        new_mixture.fractions = fractions;
        new_mixture.ratio = ratio;

        input_substrate_iterator->mixtures.push_back(new_mixture);
    }
}


void ParserOpenFlux::ParseReactions() {
    const std::string reactions_path = path_ + "/model.csv";

    std::ifstream reaction_file(reactions_path);

    // skip table head-line
    reaction_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    getline(reaction_file, raw_line);
    try {
        while (!raw_line.empty()) {
            Reaction new_reaction = FillReaction(raw_line);
            reactions_.emplace_back(new_reaction);
            if (reaction_file.eof()) {
                break;
            } else {
                getline(reaction_file, raw_line);
            }
        }
    } catch (std::runtime_error &parser_error) {
        std::cerr << "Line: " << reactions_.size() << std::endl;
        throw (parser_error);
    }

    reactions_ = modelling_utills::SortReactionsByType(reactions_);
}


std::vector<std::string> ParserOpenFlux::ReadEachLine(const std::string &path) {
    std::vector<std::string> lines;

    std::ifstream input(path);
    std::string next_line;
    while (getline(input, next_line)) {
        lines.emplace_back(next_line);
    }

    return lines;
}


std::string ParserOpenFlux::GetCell(std::stringstream &line) {
    std::string new_cell;
    getline(line, new_cell, csv_delimiter_);
    return new_cell;
}


Reaction ParserOpenFlux::FillReaction(const std::string& raw_line) {
    Reaction reaction;
    std::stringstream line(raw_line);
    reaction.id = reactions_.size();
    reaction.name = GetCell(line);
    reaction.chemical_equation = ParseChemicalEquation(line);
    GetCell(line); // we don't care about Rate
    reaction.type = ParseReactionType(GetCell(line));
    auto[new_basis, is_basis_x] = ParseBasis(GetCell(line));
    reaction.basis = new_basis;
    reaction.is_set_free = is_basis_x;
    reaction.deviation = ParseDeviation(GetCell(line));

    return reaction;
}


ChemicalEquation ParserOpenFlux::ParseChemicalEquation(std::stringstream &line) {
    const std::string substrate_equation = GetCell(line);
    const size_t substrate_delimiter_position = substrate_equation.find(reaction_side_delimiter_);

    const std::string atom_equation = GetCell(line);
    size_t atom_delimiter_position = atom_equation.find(reaction_side_delimiter_);

    bool is_atom_equation_ok = (atom_delimiter_position != std::string::npos) || atom_equation.empty();

    if (substrate_delimiter_position != std::string::npos && is_atom_equation_ok) {
        const std::string left_side_substrate_equation = substrate_equation.substr(0, substrate_delimiter_position);
        const std::string left_side_atom_equation = atom_equation.empty() ? "" :
                                                    atom_equation.substr(0, atom_delimiter_position);

        ChemicalEquationSide left = FillEquationSide(
            left_side_substrate_equation, left_side_atom_equation);

        const std::string right_side_substrate_equation = substrate_equation.substr(
            substrate_delimiter_position + 2);
        const std::string right_side_atom_equation = atom_equation.empty() ? "" :
                                                     atom_equation.substr(atom_delimiter_position + 2);

        ChemicalEquationSide right = FillEquationSide(
            right_side_substrate_equation, right_side_atom_equation);

        return ChemicalEquation{left, right};
    } else {
        throw std::runtime_error("There is reaction without = sign!");
    }

}


ChemicalEquationSide
ParserOpenFlux::FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation) {
    ChemicalEquationSide equation_side = ParseSubstrateEquationSide(substrate_equation);
    ParseAtomEquationSide(atom_equation, &equation_side);
    return equation_side;
}


ChemicalEquationSide ParserOpenFlux::ParseSubstrateEquationSide(const std::string &raw_equation) {
    ChemicalEquationSide result;
    std::stringstream equation{raw_equation};
    bool previous_token_is_coefficient{false}; // true, if previous iteration have found a coefficient
    SubstrateCoefficient last_coefficient{}; // contains previous coefficient, if previous token is a coefficient
    std::string token;
    getline(equation, token, ' ');
    while (!token.empty()) {
        try {

            // We are trying to convert token into a double
            // std::stod throws an exception when token is a substrate name
            std::size_t position_of_not_number; // see second argument of std::stod
            last_coefficient = std::stod(token, &position_of_not_number);

            if (position_of_not_number != token.size()) { // because std::stod won't throw when substrate name
                throw std::invalid_argument("");          // starts with digits
            }

            // We found a coefficient
            if (!previous_token_is_coefficient) {
                previous_token_is_coefficient = true;
            } else {
                throw std::runtime_error("There is reaction with two coefficient in a row!");
            }
        } catch (std::invalid_argument) {
            if (token != substrate_delimiter_) {
                if (!previous_token_is_coefficient) {
                    last_coefficient = 1.0;
                }
                previous_token_is_coefficient = false;
                Substrate new_substrate;
                new_substrate.name = token;
                new_substrate.coefficient = last_coefficient;
                result.push_back(new_substrate);
            }
        }
        if (equation.eof()) {
            break;
        } else {
            getline(equation, token, ' ');
        }
    }
    if (!previous_token_is_coefficient) {
        if (!result.empty()) {
            return result;
        } else {
            throw std::runtime_error("There is reaction with empty side!");
        }
    } else {
        throw std::runtime_error("There is reaction with a coefficient without substrate!");
    }
}


void ParserOpenFlux::ParseAtomEquationSide(const std::string &atom_equation, ChemicalEquationSide *equation_side) {
    // It's ok if atom_equation is empty
    if (!atom_equation.empty()) {
        std::stringstream equation{atom_equation};
        std::string token;
        getline(equation, token, ' ');
        auto current_substance = equation_side->begin();
        while (!token.empty()) {

            // We don't check coefficients
            if (token != substrate_delimiter_ && !std::isdigit(token[0])) {
                current_substance->formula = token;
                ++current_substance;
            }
            if (equation.eof()) {
                break;
            } else {
                getline(equation, token, ' ');
            }
        }
        if (current_substance != equation_side->end()) {
            throw std::runtime_error("There is contradiction between substance and atom reactions!");
        }
    }
}


ReactionType ParserOpenFlux::ParseReactionType(const std::string &type) {
    if (!type.empty()) {
        if (type == "F") {
            return ReactionType::Irreversible;
        } else if (type == "FR") {
            return ReactionType::Forward;
        } else if (type == "R") {
            return ReactionType::Backward;
        } else if (type == "S") {
            return ReactionType::IsotopomerBalance;
        } else if (type == "B") {
            return ReactionType::MetaboliteBalance;
        } else {
            throw std::runtime_error("There is reaction with bad type!");
        }
    } else {
        throw std::runtime_error("There is reaction without type in model!");
    }
}


std::tuple<Basis, bool> ParserOpenFlux::ParseBasis(const std::string &basis) {
    if (basis == "X" || basis == "x") {
        return {std::numeric_limits<double>::quiet_NaN(), true};
    } else if (basis.empty()) {
        return {std::numeric_limits<double>::quiet_NaN(), false};
    } else {
        return {std::stod(basis), true};
    }
}


Deviation ParserOpenFlux::ParseDeviation(const std::string &deviation) {
    if (deviation.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no deviation
    } else {
        return std::stod(deviation);
    }
}

} //namespace khnum
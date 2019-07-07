#include "parser_open_flux.h"

#include "sort_reactions.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>
#include <sstream>
#include <tuple>

void ParserOpenFlux::ReadExcludedMetabolites() {
    const std::string excluded_metabolites_path = path_ + "/excluded_metabolites.txt";
    std::vector<std::string> excluded_metabolites = ParseEachLine(excluded_metabolites_path);
    results_.excluded_metabolites.emplace(excluded_metabolites);
}


void ParserOpenFlux::ReadMeasuredIsotopes() {
    const std::string measured_isotopes_path = path_ + "/measured_isotopes.txt";

    std::vector<Emu> measured_isotopes;
    std::vector<std::string> raw_measured_isotopes = ParseEachLine(measured_isotopes_path);
    int line_number = 1;
    for (const std::string &raw_measured_isotope : raw_measured_isotopes) {
        std::stringstream line(raw_measured_isotope);
        Emu new_emu;
        getline(line, new_emu.name, ':');

        std::string row_atom_states;
        getline(line, row_atom_states);

        new_emu.atom_states.resize(row_atom_states.size());
        for (int atom_position = 0; atom_position < row_atom_states.size(); ++atom_position) {
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

        measured_isotopes.push_back(new_emu);
    }

    results_.measuredEmu.emplace(measured_isotopes);
}


void ParserOpenFlux::ReadMeasurements() {
    const std::string measurements_path = path_ + "/measurements.csv";

    std::vector<Measurement> measurements;
    std::ifstream input(measurements_path);
    // skip table head-line
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    for (const Emu &isotope : *results_.measuredEmu) {
        Measurement new_measurement;
        new_measurement.emu = isotope;
        for (int mass_shift = 0; mass_shift < isotope.atom_states.size() + 1; ++mass_shift) {
            getline(input, raw_line);
            std::stringstream line(raw_line);

            double value = std::stod(GetCell(line));
            new_measurement.mid.push_back(value);

            double error = std::stod(GetCell(line));
            new_measurement.errors.push_back(error);
        }
        measurements.push_back(new_measurement);
    }

    results_.measurements.emplace(measurements);
}


void ParserOpenFlux::ReadReactions() {
    const std::string reactionsPath = path_ + "/model.csv";
    std::vector<Reaction> reactions = ParseReactions(reactionsPath);
    reactions = SortReactionsByType(reactions);
    results_.reactions.emplace(reactions);
}


void ParserOpenFlux::ReadSubstrateInput() {
    const std::string input_substrates_path = path_ + "/substrate_input.csv";

    std::vector<InputSubstrate> input_substrates;
    std::ifstream input(input_substrates_path);

    // skip table head-line
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    while (getline(input, raw_line)) {
        std::stringstream line(raw_line);
        std::string input_substrate_name;
        getline(line, input_substrate_name, csv_delimiter);

        auto input_substrate_iterator = std::find_if(input_substrates.begin(),
                                                     input_substrates.end(),
                                                     [&input_substrate_name](InputSubstrate &input_substrate) {
                                                         return input_substrate.name == input_substrate_name;
                                                     });

        if (input_substrate_iterator == input_substrates.end()) {
            InputSubstrate new_input_substrate;
            new_input_substrate.name = input_substrate_name;
            input_substrates.push_back(new_input_substrate);
            input_substrate_iterator = input_substrates.end();
            --input_substrate_iterator;
        }

        // now input_substrate_iterator point to the entering substrate
        std::vector<Fraction> fractions;
        std::string raw_labeling_pattern;
        getline(line, raw_labeling_pattern, csv_delimiter);

        std::stringstream labeling_pattern(raw_labeling_pattern);
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

    results_.input_substrate.emplace(input_substrates);
}


std::vector<Reaction> ParserOpenFlux::ParseReactions(const std::string &model_path) {
    std::ifstream model_file(model_path);
    std::vector<Reaction> reactions;

    // skip table head-line
    model_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    getline(model_file, raw_line);
    int reaction_id = 0;
    try {
        while (!raw_line.empty()) {
            std::stringstream line(raw_line);
            Reaction new_reaction;
            new_reaction.id = reaction_id;
            ++reaction_id;
            FillReaction(&new_reaction, line);
            reactions.emplace_back(new_reaction);
            if (model_file.eof()) {
                break;
            } else {
                getline(model_file, raw_line);
            }
        }
    } catch (std::runtime_error &parser_error) {
        std::cerr << "Line: " << reaction_id << std::endl;
        throw (parser_error);
    }


    return reactions;
}


void ParserOpenFlux::FillReaction(Reaction *reaction, std::stringstream &line) {

    reaction->name = GetCell(line);
    reaction->chemical_equation = ParseChemicalEquation(line);
    reaction->rate = ParseRate(GetCell(line));
    reaction->type = ParseReactionType(GetCell(line));
    auto[new_basis, is_basis_x] = ParseBasis(GetCell(line));
    reaction->basis = new_basis;
    reaction->is_set_free = is_basis_x;
    reaction->deviation = ParseDeviation(GetCell(line));

    // turn this code ON if UB and LB are available
    /*
    reaction->setted_lower_bound = ParseLowerBound(GetCell(line));
    reaction->setted_upper_bound = ParseUpperBound(GetCell(line));
     */

    reaction->setted_lower_bound = 0;
    reaction->setted_upper_bound = 10;
    GetCell(line);
    GetCell(line);
}


ChemicalEquation ParserOpenFlux::ParseChemicalEquation(std::stringstream &line) {
    const std::string substrate_equation = GetCell(line);
    const size_t substrate_delimiter_position = substrate_equation.find(reaction_side_delimiter);

    const std::string atom_equation = GetCell(line);
    size_t atom_delimiter_position = atom_equation.find(reaction_side_delimiter);

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
    ParseAtomEquationSide(&equation_side, atom_equation);
    return equation_side;
}


ChemicalEquationSide ParserOpenFlux::ParseSubstrateEquationSide(const std::string &raw_equation) {
    ChemicalEquationSide result;
    std::stringstream equation{raw_equation};
    bool previous_token_is_coefficient{false}; // true, if previous iteration have found a coefficient
    SubstrateCoefficient last_coefficient{}; // contain previous coefficient, if previous token is a coefficient
    std::string token;
    getline(equation, token, ' ');
    while (!token.empty()) {
        try {
            std::size_t position_of_not_number; // see second argument of std::stod
            last_coefficient = std::stod(token, &position_of_not_number);

            if (position_of_not_number != token.size()) { // in case of substrate name starts with a number
                throw std::invalid_argument("");
            }
            if (!previous_token_is_coefficient) {
                previous_token_is_coefficient = true;

            } else {
                throw std::runtime_error("There is reaction with two coefficient in row!");
            }
        } catch (std::invalid_argument) {
            if (token != substrate_delimiter) {
                if (!previous_token_is_coefficient) {
                    last_coefficient = 1.0;
                }
                previous_token_is_coefficient = false;
                Substrate new_substrate;
                new_substrate.name = token;
                new_substrate.coefficient = last_coefficient;
                new_substrate.formula = "";
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
        throw std::runtime_error("There is reaction with coefficient without a substrate!");
    }
}


void ParserOpenFlux::ParseAtomEquationSide(ChemicalEquationSide *equation_side, const std::string &atom_equation) {
    if (!atom_equation.empty()) {
        std::stringstream equation{atom_equation};
        std::string token;
        getline(equation, token, ' ');
        auto current_substance = equation_side->begin();
        while (!token.empty()) {
            if (token != substrate_delimiter && !std::isdigit(token[0])) {
                current_substance->formula = token;
                ++current_substance;
            }
            if (equation.eof()) {
                break;
            } else {
                getline(equation, token, ' ');
            }
        }
        if (current_substance == equation_side->end()) {
            return;
        } else {
            throw std::runtime_error("There is contradiction between substance and atom reactions!");
        }
    } else { // reaction without tracer
        return;
    }
}


Rate ParserOpenFlux::ParseRate(const std::string &rate) {
    if (rate.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no rate
    } else {
        return std::stod(rate);
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


Flux ParserOpenFlux::ParseLowerBound(const std::string &lower_bound) {
    return std::stod(lower_bound);
}


Flux ParserOpenFlux::ParseUpperBound(const std::string &upper_bound) {
    return std::stod(upper_bound);
}


Deviation ParserOpenFlux::ParseDeviation(const std::string &deviation) {
    if (deviation.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no deviation
    } else {
        return std::stod(deviation);
    }
}


std::vector<std::string> ParserOpenFlux::ParseEachLine(const std::string &path) {
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
    getline(line, new_cell, csv_delimiter);
    return new_cell;
}
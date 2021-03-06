#include "parser/open_flux_parser/open_flux_utills.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <utilities/debug_utills/debug_prints.h>
#include "utilities/matrix.h"


namespace khnum {
namespace open_flux_parser {
std::vector<std::string> GetLines(const std::string &path) {
    std::vector<std::string> lines;

    std::ifstream input(path);
    if (!input) {
        throw(std::runtime_error("Can't open file: " + path));
    }

    std::string next_line;
    while (getline(input, next_line)) {
        if (!next_line.empty()) {
            lines.emplace_back(next_line);
        }
    }

    return lines;
}

std::vector<Emu> ParseMeasuredIsotopes(const std::vector<std::string> &raw_measured_isotopes) {
    std::vector<Emu> measured_isotopes;

    // line_number need for error messages
    int line_number = 1;
    for (const std::string &raw_measured_isotope : raw_measured_isotopes) {
        Emu measured_isotope = ParseOneMeasuredIsotope(raw_measured_isotope, line_number);
        measured_isotopes.emplace_back(measured_isotope);
        ++line_number;
    }

    return measured_isotopes;
}

Emu ParseOneMeasuredIsotope(const std::string& raw_measured_isotope, int line_number /* = -1 */) {
    std::stringstream line(raw_measured_isotope);
    Emu measured_isotope;
    getline(line, measured_isotope.name, ':');

    std::string row_atom_states;
    getline(line, row_atom_states);

    measured_isotope.atom_states.resize(row_atom_states.size());
    for (size_t atom_position = 0; atom_position < row_atom_states.size(); ++atom_position) {
        if (row_atom_states[atom_position] == '1') {
            measured_isotope.atom_states[atom_position] = true;
        } else if (row_atom_states[atom_position] == '0') {
            measured_isotope.atom_states[atom_position] = false;
        } else {
            throw std::runtime_error(
                "Line: " + std::to_string(line_number) + " There is strange atom states in measured isotope!");
        }
    }

    return measured_isotope;
}

std::string GetMeasuredIsotopeName(const Emu& emu) {
    std::string name = emu.name;
    name += ":";
    for (char state : emu.atom_states) {
        if (state) {
            name += "1";
        } else {
            name += "0";
        }
    }
    return name;
}

Matrix ParseCorrectionMatrix(const std::vector<std::string>& raw_matrix, const Delimiters& delimiters) {
    std::vector<double> first_line;
    std::stringstream line(raw_matrix.at(0));
    while (line) {
        try {
            double value = std::stod(GetCell(line, delimiters));
            first_line.push_back(value);
        } catch (std::invalid_argument& err) {
            break;
        }

    }
    Matrix correction_matrix(first_line.size(), first_line.size());
    for (int i = 0; i < raw_matrix.size(); ++i) {
        std::stringstream row(raw_matrix.at(i));
        int j = 0;
        while (row) {
            try {
                double value = std::stod(GetCell(row, delimiters));
                correction_matrix(i, j) = value;
                ++j;
            } catch (std::invalid_argument& err) {
                j = 0;
                break;
            }


        }
    }
    return correction_matrix;
}

std::vector<Measurement> ParseMeasurements(const std::vector<std::string>& raw_measurements,
                                           const std::vector<Emu>& measured_isotopes,
                                           const Delimiters& delimiters) {
    std::vector<Measurement> measurements;

    size_t line_number = 1;
    for (const Emu &isotope : measured_isotopes) {
        Measurement new_measurement;
        new_measurement.emu = isotope;
        for (size_t mass_shift = 0; mass_shift < isotope.atom_states.size() + 1; ++mass_shift) {
            std::stringstream line(raw_measurements.at(line_number));
            ++line_number;

            double value = std::stod(GetCell(line, delimiters));
            new_measurement.mid.push_back(value);

            double error = std::stod(GetCell(line, delimiters));
            new_measurement.errors.push_back(error);
        }
        measurements.push_back(new_measurement);
    }

    return measurements;
}

std::string GetCell(std::stringstream &line, const Delimiters& delimiters) {
    std::string new_cell;
    getline(line, new_cell, delimiters.csv_delimiter);
    return new_cell;
}


std::vector<InputSubstrate> ParseInputSubstrates(const std::vector<std::string>& raw_input_substrates,
                                                 const Delimiters& delimiters) {
    std::vector<InputSubstrate> input_substrates;

    for (size_t line_number = 1; line_number < raw_input_substrates.size(); ++line_number) {
        std::stringstream line(raw_input_substrates.at(line_number));
        std::string input_substrate_name;
        getline(line, input_substrate_name, delimiters.csv_delimiter);

        auto input_substrate_iterator = std::find_if(input_substrates.begin(),
                                                     input_substrates.end(),
                                                     [&input_substrate_name](InputSubstrate &input_substrate) {
                                                         return input_substrate.name == input_substrate_name;
                                                     });

        if (input_substrate_iterator == input_substrates.end()) {
            InputSubstrate new_input_substrate;
            new_input_substrate.name = input_substrate_name;
            input_substrates.push_back(new_input_substrate);

            //move iterator back, so we can work with the new substrate
            input_substrate_iterator = input_substrates.end();
            --input_substrate_iterator;
        }

        std::string raw_labeling_pattern;
        getline(line, raw_labeling_pattern, delimiters.csv_delimiter);
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

    return input_substrates;
}




std::vector<Reaction> ParseReactions(const std::string& reactions_path,
                                     const Delimiters& delimiters) {
    std::vector<Reaction> reactions;
    std::ifstream reaction_file(reactions_path);

    // skip table head-line
    reaction_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    getline(reaction_file, raw_line);
    int current_id = 0;
    while (!raw_line.empty()) {
        Reaction new_reaction = FillReaction(raw_line, current_id, delimiters);
        ++current_id;
        reactions.emplace_back(new_reaction);
        if (reaction_file.eof()) {
            break;
        } else {
            getline(reaction_file, raw_line);
        }
    }
    return reactions;
}

Reaction FillReaction(const std::string& raw_line, int id, const Delimiters& delimiters) {
    Reaction reaction;
    std::stringstream line(raw_line);
    reaction.id = id;
    reaction.name = GetCell(line, delimiters);
    reaction.chemical_equation = ParseChemicalEquation(line, delimiters);
    GetCell(line, delimiters); // we don't care about Rate
    reaction.type = ParseReactionType(GetCell(line, delimiters));
    std::optional<Basis> basis = ParseBasis(GetCell(line, delimiters));
    if (basis) {
        reaction.basis = *basis;
        reaction.is_set_free = true;
    } else {
        reaction.basis = std::numeric_limits<double>::quiet_NaN();
        reaction.is_set_free = false;
    }
    reaction.deviation = ParseDeviation(GetCell(line, delimiters));

    return reaction;
}


ChemicalEquation ParseChemicalEquation(std::stringstream &line, const Delimiters& delimiters) {
    const std::string substrate_equation = GetCell(line, delimiters);
    const size_t reaction_delimiter_position = substrate_equation.find(delimiters.reaction_side_delimiter);

    const std::string atom_equation = GetCell(line, delimiters);
    size_t atom_delimiter_position = atom_equation.find(delimiters.reaction_side_delimiter);

    bool is_atom_equation_ok = (atom_delimiter_position != std::string::npos) || atom_equation.empty();

    if (reaction_delimiter_position != std::string::npos && is_atom_equation_ok) {
        ChemicalEquation result;
        const std::string left_side_substrate_equation = substrate_equation.substr(0, reaction_delimiter_position);
        const std::string left_side_atom_equation = atom_equation.empty() ? "" :
                                                    atom_equation.substr(0, atom_delimiter_position);

        result.left = FillEquationSide(
            left_side_substrate_equation, left_side_atom_equation, delimiters);

        const std::string right_side_substrate_equation = substrate_equation.substr(
            reaction_delimiter_position + 2);
        const std::string right_side_atom_equation = atom_equation.empty() ? "" :
                                                     atom_equation.substr(atom_delimiter_position + 2);

        result.right = FillEquationSide(
            right_side_substrate_equation, right_side_atom_equation, delimiters);

        result.atom_transitions = ParseAtomTransitions(atom_equation, delimiters, &result.left, &result.right);

        return result;
    } else {
        throw std::runtime_error("There is reaction without = sign!");
    }

}


ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation,
                                      const Delimiters& delimiters) {
    ChemicalEquationSide equation_side = ParseSubstrateEquationSide(substrate_equation, delimiters);
    return equation_side;
}


ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation, const Delimiters& delimiters) {
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
            double coefficient = std::stod(token, &position_of_not_number);

            if (position_of_not_number != token.size()) { // because std::stod won't throw when substrate name
                throw std::invalid_argument("");          // starts with digits
            }

            last_coefficient = coefficient;

            // We found a coefficient
            if (!previous_token_is_coefficient) {
                previous_token_is_coefficient = true;
            } else {
                throw std::runtime_error("There is reaction with two coefficient in a row!");
            }
        } catch (std::invalid_argument&) {
            if (token != delimiters.substrate_delimiter) {
                if (!previous_token_is_coefficient) {
                    last_coefficient = 1.0;
                }
                previous_token_is_coefficient = false;
                Substrate new_substrate;
                new_substrate.id = result.size();
                new_substrate.name = token;
                new_substrate.substrate_coefficient_ = last_coefficient;
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


std::vector<AtomTransition> ParseAtomTransitions(const std::string &atom_equation, const Delimiters& delimiters, ChemicalEquationSide *left_side, ChemicalEquationSide *right_side) {
    std::vector<AtomTransition> result;

    // It's ok if atom_equation is empty
    if (atom_equation.empty()) {
        return result;
    }

    std::vector<std::string> left;
    std::vector<std::string> right;
    bool is_left = true;

    std::stringstream equation{atom_equation};
    std::string token;
    bool previous_token_is_coefficient{false}; // true, if previous iteration have found a coefficient
    SubstrateCoefficient last_coefficient{}; // contains previous coefficient, if previous token is a coefficient

    auto current_substance = left_side->begin();

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
        } catch (std::invalid_argument&) {
            if (token != delimiters.substrate_delimiter) {
                if (token != delimiters.reaction_side_delimiter) {
                    if (!previous_token_is_coefficient) {
                        last_coefficient = 1.0;
                    }
                    previous_token_is_coefficient = false;
                    current_substance->size = token.size();
                    ++current_substance;
                    if (is_left) {
                        left.push_back(token);
                    } else {
                        right.push_back(token);
                    }
                } else {
                    is_left = false;
                    previous_token_is_coefficient = false;
                    current_substance = right_side->begin();
                }
            }
        }
        if (equation.eof()) {
            break;
        } else {
            getline(equation, token, ' ');
        }
    }

    for (int substrate_pos = 0; substrate_pos < left.size(); ++substrate_pos) {
        const std::string &substrate = left.at(substrate_pos);
        for (int substrate_atom = 0; substrate_atom < substrate.size(); ++substrate_atom) {
            AtomTransition transition;
            transition.substrate_pos = substrate_pos;
            transition.substrate_atom = substrate_atom;

            for (int product_pos = 0; product_pos < right.size(); ++product_pos) {
                const std::string &product = right.at(product_pos);
                for (int product_atom = 0; product_atom < product.size(); ++product_atom) {
                    if (substrate[substrate_atom] == product[product_atom]) {
                        transition.product_pos = product_pos;
                        transition.product_atom = product_atom;
                        result.push_back(transition);
                        break;
                    }
                }
            }


        }
    }

    return result;
}


ReactionType ParseReactionType(const std::string &type) {
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


std::optional<Basis> ParseBasis(const std::string &basis) {
    if (basis == "X" || basis == "x") {
        return {std::numeric_limits<double>::quiet_NaN()};
    } else if (basis.empty()) {
        return {};
    } else {
        return {std::stod(basis)};
    }
}


Deviation ParseDeviation(const std::string &deviation) {
    if (deviation.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no deviation
    } else {
        return std::stod(deviation);
    }
}

} // namespace open_flux_parser
} // namespace khnum
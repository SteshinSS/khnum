#include "model_parser.h"
#include "parser_utilities.h"
#include <exception>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>

std::vector<Reaction> ParseReactions(const std::string &model_path) {
    std::ifstream model_file(model_path);
    std::vector<Reaction> reactions;

    // skip table head-line
    model_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    getline(model_file, raw_line);
    while (!raw_line.empty()) {
        std::stringstream line(raw_line);
        Reaction new_reaction;
        FillReaction(&new_reaction, line);
        reactions.emplace_back(new_reaction);
        if (model_file.eof()) {
            break;
        } else {
            getline(model_file, raw_line);
        }
    }

    return reactions;
}

void FillReaction(Reaction *reaction, std::stringstream &line) {
    reaction->name = GetCell(line);
    reaction->chemical_equation = ParseChemicalEquation(line);
    reaction->rate = ParseRate(GetCell(line));
    reaction->type = ParseReactionType(GetCell(line));
    auto [new_basis, is_basis_x] = ParseBasis(GetCell(line));
    reaction->basis = new_basis;
    reaction->is_basis_x = is_basis_x;
    reaction->deviation = ParseDeviation(GetCell(line));
    reaction->lower_bound = ParseLowerBound(GetCell(line));
    reaction->upper_bound = ParseUpperBound(GetCell(line));
    reaction->is_prime_basis = ParsePrimeBasis(line);
}

ChemicalEquation ParseChemicalEquation(std::stringstream &line) {
    const std::string substrate_equation = GetCell(line);
    const std::string atom_equation = GetCell(line);
    size_t substrate_delimiter_position = substrate_equation.find(reaction_side_delimiter);
    size_t atom_delimiter_position = atom_equation.find(reaction_side_delimiter);
    bool is_atom_equation_ok = (atom_delimiter_position != std::string::npos) || atom_equation.empty();

    if (substrate_delimiter_position != std::string::npos && is_atom_equation_ok) {
        const std::string left_side_substrate_equation = substrate_equation.substr(0, substrate_delimiter_position);
        const std::string left_side_atom_equation = atom_equation.substr(0, atom_delimiter_position);
        ChemicalEquationSide left = FillEquationSide(
            left_side_substrate_equation, left_side_atom_equation);

        const std::string right_side_substrate_equation = substrate_equation.substr(substrate_delimiter_position + 2);
        const std::string right_side_atom_equation = atom_equation.empty() ? "" :
                                                     atom_equation.substr(atom_delimiter_position + 2);

        ChemicalEquationSide right = FillEquationSide(
            right_side_substrate_equation, right_side_atom_equation);

        return ChemicalEquation{left, right};
    } else {
        throw std::runtime_error("There is reaction without = sign!");
    }

}

ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation) {
    ChemicalEquationSide equation_side = ParseSubstrateEquationSide(substrate_equation);
    ParseAtomEquationSide(&equation_side, atom_equation);
    return equation_side;
}

ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation) {
    ChemicalEquationSide result;
    std::stringstream equation{raw_equation};
    bool previous_token_is_coefficient{false}; // true, if previous iteration have found a coefficient
    SubstrateCoefficient last_coefficient{}; // contain previous coefficient, if previous token is coefficient
    std::string token;
    getline(equation, token, ' ');
    while (!token.empty()) {
        if (std::isdigit(token[0])) { // it is coefficient
            if (!previous_token_is_coefficient) {
                previous_token_is_coefficient = true;
                last_coefficient = std::stod(token);
            } else {
                throw std::runtime_error("There is reaction with two coefficient in row!");
            }
        } else if (token != substrate_delimiter) {
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

void ParseAtomEquationSide(ChemicalEquationSide *equation_side, const std::string &atom_equation) {
    if (!atom_equation.empty()) {
        std::stringstream equation{atom_equation};
        std::string token;
        getline(equation, token, ' ');
        ChemicalEquationSide::iterator current_substance = equation_side->begin();
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

Rate ParseRate(const std::string &rate) {
    if (rate.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no rates
    } else {
        return std::stod(rate);
    }
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
            throw std::runtime_error("There is reaction with bad Type!");
        }
    } else {
        throw std::runtime_error("There is reaction without Type in model!");
    }
}

std::tuple<Basis, bool> ParseBasis(const std::string &basis) {
    if (basis == "X" || basis == "x") {
        return {std::numeric_limits<double>::quiet_NaN(), true};
    } else if (basis.empty()) {
        return {std::numeric_limits<double>::quiet_NaN(), false};
    }
    else {
        return {std::stod(basis), false};
    }
}

Bound ParseLowerBound(const std::string &lower_bound) {
    return std::stod(lower_bound);
}

Bound ParseUpperBound(const std::string &upper_bound) {
    return std::stod(upper_bound);
}

Deviation ParseDeviation(const std::string &deviation) {
    if (deviation.empty()) {
        return std::numeric_limits<double>::quiet_NaN();  // no deviation
    } else {
        return std::stod(deviation);
    }
}

bool ParsePrimeBasis(std::stringstream &line) {
    std::string last_cell;
    getline(line, last_cell, '\r');
    if (last_cell.empty()) {
        return false;
    } else {
        if (last_cell == "*") { // warning
            return true;
        } else {
            throw std::runtime_error(
                "There is reaction with strange last cell. Is asterisk(*) ok?");
        }
    }
}

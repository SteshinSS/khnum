#include "parser/python_parser.h"

#include <iostream>
#include <sstream>

#include "utilities/reaction.h"

namespace khnum {
ParserResults ParserMaranas::GetResults() {
    ParserResults results;
    results.reactions = reactions_;
    results.measured_isotopes = measured_isotopes_;
    results.measurements = measurements_;
    results.input_substrate = input_substrates_;
    results.excluded_metabolites = excluded_metabolites_;

    return results;
}

void ParserMaranas::ParseReactions() {
    Py_Initialize();
    PyRun_SimpleString("import os");
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(os.getcwd() + \"/../src/parser/maranas_parser\")");

    PyObject *pName = PyUnicode_DecodeFSDefault("Maranas");
    PyObject *pModule = PyImport_Import(pName);
    if (!pModule) {
        exit(1);
    }
    PyObject *pFunc = PyObject_GetAttrString(pModule, "parse");
    if (!pFunc) {
        exit(2);
    }
    PyObject *pValue = PyObject_CallObject(pFunc, NULL);
    if (!pValue) {
        exit(3);
    }
    const char *command = PyBytes_AsString(pValue);

    std::string answer(command);
    std::stringstream stream(answer);
    int total_reactions = 0;
    stream >> total_reactions;
    std::cout << total_reactions << std::endl;
    std::cout << answer << std::endl;
    for (int i = 0; i < total_reactions; ++i) {
        reactions_.push_back(ParseReaction(stream));
    }
    for (int i = 0; i < total_reactions; ++i) {
        const Reaction &reaction = reactions_[i];
        if (reaction.type == ReactionType::Forward) {
            Reaction reversed(reaction);
            reversed.type = ReactionType::Backward;
            reversed.id = reactions_.size();
            reactions_.push_back(reversed);
        }
    }
}

Reaction ParserMaranas::ParseReaction(std::stringstream &stream) {
    Reaction result;
    stream >> result.id;
    stream >> result.name;
    bool is_reversed = false;
    stream >> is_reversed;
    bool is_excluded = false;
    stream >> is_excluded;

    ChemicalEquationSide left;
    int total_left = 0;
    stream >> total_left;
    for (int i = 0; i < total_left; ++i) {
        Substrate sub;
        stream >> sub.id;
        stream >> sub.substrate_coefficient_;
        stream >> sub.name;
        stream >> sub.size;
        left.push_back(sub);
    }

    ChemicalEquationSide right;
    int total_right = 0;
    stream >> total_right;
    for (int i = 0; i < total_right; ++i) {
        Substrate sub;
        stream >> sub.id;
        stream >> sub.substrate_coefficient_;
        stream >> sub.name;
        stream >> sub.size;
        right.push_back(sub);
    }

    result.chemical_equation.left = left;
    result.chemical_equation.right = right;
    int total_transitions = 0;
    stream >> total_transitions;
    for (int i = 0; i < total_transitions; ++i) {
        AtomTransition transition;
        stream >> transition.substrate_pos >> transition.substrate_atom;
        stream >> transition.product_pos >> transition.product_atom;
        result.chemical_equation.atom_transitions.push_back(transition);
    }

    if (is_excluded) {
        result.type = ReactionType::MetaboliteBalance;
    } else {
        if (!is_reversed) {
            result.type = ReactionType::Irreversible;
        } else {
            result.type = ReactionType::Forward;
        }
    }

    result.basis = std::numeric_limits<double>::quiet_NaN();
    result.is_set_free = false;
    result.deviation = std::numeric_limits<double>::quiet_NaN();
    return result;
}

void ParserMaranas::ParseExcludedMetabolites() {

}

void ParserMaranas::ParseMeasuredIsotopes() {

}

void ParserMaranas::ParseMeasurements() {

}

void ParserMaranas::ParseCorrectionMatrices() {

}

void ParserMaranas::ParseSubstrateInput() {

}

}


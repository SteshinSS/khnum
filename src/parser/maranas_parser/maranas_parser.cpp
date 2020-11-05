#include "parser/maranas_parser.h"

#include <iostream>
#include <sstream>

#include "utilities/reaction.h"
#include "parser/open_flux_parser/open_flux_utills.h"
#include "modeller/sort_reactions.h"
#include "utilities/debug_utills/debug_prints.h"

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

void ParserMaranas::Parse() {
    ParseReactions();
    ParseExcludedMetabolites();
    ParseMeasurements();
    ParseSubstrateInput();
    PrintReactions(reactions_);
    reactions_ = modelling_utills::SortReactionsByType(reactions_);

}

void ParserMaranas::ParseReactions() {
    Py_Initialize();

    // change current directory
    PyRun_SimpleString("import os");
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(os.getcwd() + \"/../src/parser/maranas_parser\")");

    PyObject *pName = PyUnicode_DecodeFSDefault("Maranas");
    PyObject *pModule = PyImport_Import(pName);
    PyObject *pFunc = PyObject_GetAttrString(pModule, "parse");
    PyObject *pArgs = PyTuple_New(1);
    PyObject *path = Py_BuildValue("s", path_.c_str());
    PyTuple_SetItem(pArgs, 0, path);
    PyObject *pValue = PyObject_CallObject(pFunc, pArgs);
    const char *command = PyBytes_AsString(pValue);
    std::string answer(command);

    std::stringstream stream(answer);
    int total_reactions = 0;
    stream >> total_reactions;
    for (int i = 0; i < total_reactions; ++i) {
        reactions_.push_back(ParseReaction(stream));
    }

    // Add backward reactions for the reversible
    for (int i = 0; i < total_reactions; ++i) {
        const Reaction &reaction = reactions_[i];

        CollectSubstrateSizes(reaction.chemical_equation.left);
        CollectSubstrateSizes(reaction.chemical_equation.right);

        if (reaction.type == ReactionType::Forward) {
            Reaction reversed(reaction);
            reversed.type = ReactionType::Backward;
            reversed.id = reactions_.size();
            ChemicalEquationSide left = reversed.chemical_equation.left;
            reversed.chemical_equation.left = reversed.chemical_equation.right;
            reversed.chemical_equation.right = left;
            for (AtomTransition &transition : reversed.chemical_equation.atom_transitions) {
                int substrate_position = transition.substrate_pos;
                int substrate_atom = transition.substrate_atom;

                transition.substrate_pos = transition.product_pos;
                transition.substrate_atom = transition.product_atom;

                transition.product_pos = substrate_position;
                transition.product_atom = substrate_atom;
            }
            reactions_.push_back(reversed);
        }
    }

}

void ParserMaranas::CollectSubstrateSizes(const ChemicalEquationSide &side) {
    for (const Substrate &substrate : side) {
        if (substrate_sizes_.find(substrate.name) != substrate_sizes_.end()) {
            if (substrate.size == 0) {
                continue;
            }

            if (substrate_sizes_[substrate.name] == 0) {
                substrate_sizes_[substrate.name] = substrate.size;
            } else {
                if (substrate_sizes_[substrate.name] != substrate.size) {
                    throw std::runtime_error("There are contradictions in substrate " + substrate.name + " size.");
                }
            }

        } else {
            substrate_sizes_[substrate.name] = substrate.size;
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
    for (const Reaction& reaction : reactions_) {
        if (reaction.type == ReactionType::MetaboliteBalance) {
            excluded_metabolites_.push_back(reaction.chemical_equation.left[0].name);
        }
    }
}


void ParserMaranas::ParseMeasurements() {
    const std::string measurements_path = "../modelMaranas/measurements.csv";
    std::vector<std::string> raw_measurements = open_flux_parser::GetLines(measurements_path);
    for (int i = 1; i < raw_measurements.size(); ++i) {
        Measurement measurement;
        std::stringstream line(raw_measurements[i]);
        std::string formula;
        getline(line, formula, '-');
        measurement.emu.name = formula.substr(1) + "[d]"; // skip first " symbols

        measurement.emu.atom_states = AtomStates(substrate_sizes_[measurement.emu.name], 0);
        measurement.mid = Mid(substrate_sizes_[measurement.emu.name], 0);
        std::string atoms;
        getline(line, atoms, '"');
        std::stringstream atoms_stream(atoms);
        std::string atom;
        getline(atoms_stream, atom, ',');
        while (!atom.empty()) {
            int atom_pos = std::stoi(atom);
            measurement.emu.atom_states[atom_pos - 1] = 1;
            if (atoms_stream.eof()) {
                break;
            }
            getline(atoms_stream, atom, ',');
        }
        line.ignore(1); // skip ,
        for (int i = 0; i < measurement.mid.size(); ++i) {
            double mass = 0.0;
            line >> mass;
            measurement.mid[i] = mass;
        }
        measurement.errors = Errors(substrate_sizes_[measurement.emu.name], 0.005);
        measurements_.push_back(measurement);
        measured_isotopes_.push_back(measurement.emu);
    }
}


void ParserMaranas::ParseSubstrateInput() {

    const std::string input_substrates_path = "../modelMaranas/substrate_input.csv";
    const std::vector<std::string>& raw_substrates = open_flux_parser::GetLines(input_substrates_path);
    std::vector<InputSubstrate> input_substrates;

    for (size_t line_number = 1; line_number < raw_substrates.size(); ++line_number) {
        std::stringstream line(raw_substrates.at(line_number));
        std::string input_substrate_name;
        getline(line, input_substrate_name, ',');

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
        getline(line, raw_labeling_pattern, ',');
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

    input_substrates_ = input_substrates;
}

}


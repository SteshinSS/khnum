#include "model_supplementary_parser.h"
#include "parser_utilities.h"
#include "../utilities/measurements_struct.h"
#include "../utilities/EMU.h"
#include "../utilities/MID.h"
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>

std::vector<std::string> ParseExcludedMetabolites(const std::string &excluded_metabolites_path) {
    return ParseEachLine(excluded_metabolites_path);
}

std::vector<EMU> ParseMeasuredIsotopes(const std::string &measured_isotopes_path) {
    std::vector<EMU> measured_isotopes;
    std::vector<std::string> raw_measured_isotopes = ParseEachLine(measured_isotopes_path);
    for (std::string const &raw_measured_isotope : raw_measured_isotopes) {
        std::stringstream line(raw_measured_isotope);
        EMU new_emu;
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
                throw std::runtime_error("There is strange atom states in measured isotope!");
            }
        }

        measured_isotopes.push_back(new_emu);
    }
    return measured_isotopes;
}

std::vector<Measurement> ParseMeasurments(const std::string &measurements_path,
                                        const std::vector<EMU> &measured_isotopes) {
    std::vector<Measurement> measurements;
    std::ifstream input(measurements_path);
    // skip table head-line
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string raw_line;
    for (const EMU &isotope : measured_isotopes) {
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

    return measurements;
}

std::vector<InputSubstrate> ParseInputSubstrates(const std::string &input_substrates_path) {
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

    return input_substrates;
}

std::vector<std::string> ParseEachLine(const std::string &file_path) {
    std::vector<std::string> lines;

    std::ifstream input(file_path);
    std::string next_line;
    while (getline(input, next_line)) {
        lines.emplace_back(next_line);
    }

    return lines;
}


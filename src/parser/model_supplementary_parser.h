#ifndef CFLEX_MODEL_SUPPLEMENTARY_PARSER_H
#define CFLEX_MODEL_SUPPLEMENTARY_PARSER_H

#include "measurements_struct.h"
#include "EMU.h"
#include "input_substrate.h"
#include "MID.h"
#include <vector>
#include <string>

std::vector<std::string> ParseExcludedMetabolites(const std::string &excluded_metabolites_path);

std::vector<EMU> ParseMeasuredIsotopes(const std::string &measured_isotopes_path);

std::vector<Measurement> ParseMeasurments(const std::string &measurements_path,
                                        const std::vector<EMU> &measured_isotopes);

std::vector<InputSubstrate> ParseInputSubstrates(const std::string &input_substrates_path);

std::vector<std::string> ParseEachLine(const std::string &file_path);

#endif //CFLEX_MODEL_SUPPLEMENTARY_PARSER_H

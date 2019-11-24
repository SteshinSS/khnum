#pragma once

#include "parser/parser.h"
#include "parser/open_flux_parser/open_flux_utills.h"
#include "utilities/matrix.h"


namespace khnum {

// Parse OpenFlux[2] input files
class ParserOpenFlux : public IParser {
public:
    ParserOpenFlux(std::string path) : path_(std::move(path)) {
        delimiters_.csv_delimiter = ',';
        delimiters_.substrate_delimiter = "+";
        delimiters_.reaction_side_delimiter = "=";
    }

    ParserResults GetResults() override;

    void ParseExcludedMetabolites() override;

    void ParseMeasuredIsotopes() override;

    void ParseMeasurements() override;

    void ParseCorrectionMatrices() override;

    void ParseSubstrateInput() override;

    void ParseReactions() override;

    inline void SetCsvDelimeter(char delimiter) {
       delimiters_.csv_delimiter = delimiter;
    }

    inline void SetReactionSideDelimeter(const std::string &delimiter) {
        delimiters_.reaction_side_delimiter = delimiter;
    }

    inline void SetSubstrateDelimiter(const std::string &delimiter) {
        delimiters_.substrate_delimiter = delimiter;
    }


private:
    const std::string path_;

    std::vector<Reaction> reactions_;
    std::vector<Emu> measured_isotopes_;
    std::vector<Measurement> measurements_;
    std::vector<InputSubstrate> input_substrates_;
    std::vector<std::string> excluded_metabolites_;

    open_flux_parser::Delimiters delimiters_;
};

} //namespace khnum
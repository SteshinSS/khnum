#pragma once

#include "parser/parser.h"


namespace khnum {

// Parse OpenFlux[2] input files
class ParserOpenFlux : public Parser {
public:
    ParserOpenFlux(const std::string &path) : path_(path) {}

    ParserResults GetResults() override;

    void ParseExcludedMetabolites() override;

    void ParseMeasuredIsotopes() override;

    void ParseMeasurements() override;

    void ParseSubstrateInput() override;

    void ParseReactions() override;

    inline void SetCsvDelimeter(char delimiter) {
        csv_delimiter_ = delimiter;
    }

    inline void SetReactionSideDelimeter(const std::string &delimiter) {
        reaction_side_delimiter_ = delimiter;
    }

    inline void SetSubstrateDelimiter(const std::string &delimiter) {
        substrate_delimiter_ = delimiter;
    }


private:
    const std::string path_;

    std::vector<Reaction> reactions_;
    std::vector<Emu> measured_isotopes_;
    std::vector<Measurement> measurements_;
    std::vector<InputSubstrate> input_substrates_;
    std::vector<std::string> excluded_metabolites_;

    char csv_delimiter_{','};

    // it is string in case of "->" delimiter
    std::string reaction_side_delimiter_{"="};
    std::string substrate_delimiter_{"+"};

private:
    std::vector<std::string> ReadEachLine(const std::string &path);

    inline std::string GetCell(std::stringstream &line);

    Reaction FillReaction(const std::string& raw_line);

    ChemicalEquation ParseChemicalEquation(std::stringstream &line);

    ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation);

    ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation);

    void ParseAtomEquationSide(const std::string &atom_equation, ChemicalEquationSide *equation_side);

    Rate ParseRate(const std::string &rate);

    ReactionType ParseReactionType(const std::string &type);

    std::tuple<Basis, bool> ParseBasis(const std::string &basis);

    Deviation ParseDeviation(const std::string &deviation);

    Flux ParseLowerBound(const std::string &lower_bound);

    Flux ParseUpperBound(const std::string &upper_bound);



};

} //namespace khnum
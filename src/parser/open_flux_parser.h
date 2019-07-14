#pragma once

#include "parser/parser.h"

class ParserOpenFlux : public Parser {
public:

    ParserOpenFlux(const std::string& path) : path_(path) {}

    ParserResults GetResults() override;

    void ReadExcludedMetabolites() override;
    void ReadMeasuredIsotopes() override;
    void ReadMeasurements() override;
    void ReadReactions() override;
    void ReadSubstrateInput() override;


private:
    const std::string path_;

    std::vector<Reaction> reactions_;
    std::vector<Emu> measured_isotopes_;
    std::vector<Measurement> measurements_;
    std::vector<InputSubstrate> input_substrates_;
    std::vector<std::string> excluded_metabolites_;

    char csv_delimiter{','};
    std::string reaction_side_delimiter{"="};
    std::string substrate_delimiter{"+"};


private:
    std::vector<Reaction> ParseReactions(const std::string &model_path);

    void FillReaction(Reaction *reaction, std::stringstream &line);
    ChemicalEquation ParseChemicalEquation(std::stringstream &line);
    Rate ParseRate(const std::string &rate);
    ReactionType ParseReactionType(const std::string &type);
    std::tuple<Basis, bool> ParseBasis(const std::string &basis);
    Deviation ParseDeviation(const std::string &deviation);
    Flux ParseLowerBound(const std::string &lower_bound);
    Flux ParseUpperBound(const std::string &upper_bound);
    ChemicalEquationSide FillEquationSide(const std::string &substrate_equation, const std::string &atom_equation);
    ChemicalEquationSide ParseSubstrateEquationSide(const std::string &raw_equation);
    void ParseAtomEquationSide(ChemicalEquationSide *equation_side, const std::string &atom_equation);
    std::vector<std::string> ParseEachLine(const std::string &path);
    inline std::string GetCell(std::stringstream &line);
};

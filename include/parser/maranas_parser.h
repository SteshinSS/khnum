#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <unordered_map>
#include "parser/parser.h"


namespace khnum {

class ParserMaranas : public IParser {
 public:
  ParserMaranas(const std::string& path) : path_{path} {

  }

  void Parse() override;

  ParserResults GetResults() override;
 private:
  void ParseExcludedMetabolites();
  void ParseMeasurements();
  void ParseSubstrateInput();
  void ParseReactions();
  void CollectSubstrateSizes(const ChemicalEquationSide &side);
  bool AreThereComplexConvolutions();

  const std::string path_;
  std::unordered_map<std::string, int> substrate_sizes_;
  Reaction ParseReaction(std::stringstream &stream);

  std::vector<Reaction> reactions_;
  std::vector<Emu> measured_isotopes_;
  std::vector<Measurement> measurements_;
  std::vector<InputSubstrate> input_substrates_;
  std::vector<std::string> excluded_metabolites_;
};


} //namespace khnum
#pragma once

#include <string>

#include "parser/parser_results.h"


namespace khnum {

// Basic interface for parsers implementations
class IParser {
public:
    virtual ParserResults GetResults() = 0;

    virtual void ParseExcludedMetabolites() = 0;

    virtual void ParseMeasuredIsotopes() = 0;

    virtual void ParseMeasurements() = 0;

    virtual void ParseReactions() = 0;

    virtual void ParseSubstrateInput() = 0;

    virtual ~IParser() {};

};
} //namespace khnum
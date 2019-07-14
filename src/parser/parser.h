#pragma once

#include <string>

#include "parser/parser_results.h"


class Parser {
public:
    virtual ParserResults GetResults() = 0;

    virtual void ReadExcludedMetabolites() = 0;

    virtual void ReadMeasuredIsotopes() = 0;

    virtual void ReadMeasurements() = 0;

    virtual void ReadReactions() = 0;

    virtual void ReadSubstrateInput() = 0;

    virtual ~Parser() {};

};

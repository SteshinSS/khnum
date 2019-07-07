#ifndef CFLEX_PARSER_H
#define CFLEX_PARSER_H


#include "parser_results.h"

#include <string>

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

#endif //CFLEX_PARSER_H

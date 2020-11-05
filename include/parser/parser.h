#pragma once

#include <string>

#include "parser/parser_results.h"


namespace khnum {

// Basic interface for parsers implementations
class IParser {
public:
    virtual void Parse() = 0;

    virtual ParserResults GetResults() = 0;

    virtual ~IParser() {};
};
} //namespace khnum
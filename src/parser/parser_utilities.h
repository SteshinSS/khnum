#ifndef CFLEX_PARSER_PREFERENCES_H
#define CFLEX_PARSER_PREFERENCES_H

#include <string>
#include <sstream>

constexpr char csv_delimiter{','};
const std::string reaction_side_delimiter{"="};
const std::string substrate_delimiter{"+"};


// presume line is a whole line of csv table
// returns next cell of this table
inline std::string GetCell(std::stringstream &line) {
    std::string new_cell;
    getline(line, new_cell, csv_delimiter);
    return new_cell;
}

#endif //CFLEX_PARSER_PREFERENCES_H

#include "catch/catch.hpp"

#include <stdio.h>
#include <iostream>

#include "parser/open_flux_parser/open_flux_utills.h"

using namespace khnum;
using namespace khnum::open_flux_parser;

const std::string folder = "../tests/open_flux_parser_test";

TEST_CASE("GetLines() with empty file") {
    std::vector<std::string> result = GetLines(folder + "/empty_file.txt");
    std::vector<std::string> empty_vector = {};
    REQUIRE(result == empty_vector);
}

TEST_CASE("GetLines() with two lines") {
    std::vector<std::string> result = GetLines(folder + "/file_with_two_lines.txt");
    std::vector<std::string> lines = {"Hello", "World!"};
    REQUIRE(result == lines);
}

TEST_CASE("GetLines() with next line symbol") {
    std::vector<std::string> result = GetLines(folder + "/file_with_newline.txt");
    std::vector<std::string> lines = {"Hello", "World!"};
    REQUIRE(result == lines);
}

TEST_CASE("GetLines() with empty string") {
    std::vector<std::string> result = GetLines(folder + "/file_with_empty_line.txt");
    std::vector<std::string> lines = {"Hello", "World!"};
    REQUIRE(result == lines);
}

TEST_CASE("GetLines() with empty symbols") {
    std::vector<std::string> result = GetLines(folder + "/file_with_empty_symbols.txt");
    std::vector<std::string> line = {"Hello   world hi"};
    REQUIRE(result == line);
}
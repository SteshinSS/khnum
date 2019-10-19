#include "catch/catch.hpp"

#include <stdio.h>
#include <iostream>
#include <sstream>

#include "parser/open_flux_parser/open_flux_utills.h"

using namespace khnum;
using namespace khnum::open_flux_parser;

const std::string folder = "../tests/open_flux_parser_test";


TEST_CASE("ParseReactionType()") {
    SECTION("F") {
        std::string type = "F";
        ReactionType result = ParseReactionType(type);
        REQUIRE(result == ReactionType::Irreversible);
    }

    SECTION("FR") {
        std::string type = "FR";
        ReactionType result = ParseReactionType(type);
        REQUIRE(result == ReactionType::Forward);
    }

    SECTION("R") {
        std::string type = "R";
        ReactionType result = ParseReactionType(type);
        REQUIRE(result == ReactionType::Backward);
    }

    SECTION("S") {
        std::string type = "S";
        ReactionType result = ParseReactionType(type);
        REQUIRE(result == ReactionType::IsotopomerBalance);
    }

    SECTION("B") {
        std::string type = "B";
        ReactionType result = ParseReactionType(type);
        REQUIRE(result == ReactionType::MetaboliteBalance);
    }

    SECTION("Wrong") {
        std::string type = "FRance";
        REQUIRE_THROWS_WITH(ParseReactionType(type), "There is reaction with bad type!");
    }

    SECTION("Empty") {
        std::string type = "";
        REQUIRE_THROWS_WITH(ParseReactionType(type), "There is reaction without type in model!");
    }
}

TEST_CASE("ParseBasis()") {
    SECTION("Normal") {
        std::string basis = "1.337";
        auto result = ParseBasis(basis);
        REQUIRE(*result == Approx(1.337));
    }

    SECTION("Long") {
        std::string basis = "2.2813371337133713371337";
        auto result = ParseBasis(basis);
        REQUIRE(*result == Approx(2.2813371337));
    }

    SECTION("X") {
        std::string basis = "X";
        auto result = ParseBasis(basis);
        REQUIRE(std::isnan(*result));
    }

    SECTION("Wrong") {
        std::string basis = "XXX";
        REQUIRE_THROWS_WITH(ParseBasis(basis), "");
    }

    SECTION("Empty") {
        std::string basis = "";
        REQUIRE_THROWS_WITH(ParseBasis(basis), "");
    }
}

TEST_CASE("ParseDeviation()") {
    SECTION("Normal") {
        std::string deviation = "2.228";
        auto result = ParseDeviation(deviation);
        REQUIRE(result == Approx(2.228));
    }

    SECTION("Wrong") {
        std::string deviation = "2.aaa";
        REQUIRE_THROWS_WITH(ParseDeviation(deviation), "");
    }

    SECTION("Empty") {
        std::string deviation;
        REQUIRE_THROWS_WITH(ParseDeviation(deviation), "");
    }
}

TEST_CASE("ParseSubstrateEquationSide()") {
    Delimiters delimiters;
    delimiters.csv_delimiter = ',';
    delimiters.reaction_side_delimiter = "=";
    delimiters.substrate_delimiter = '+';

    SECTION("Long Boi") {
        std::string raw_reaction = "0.054 trp_L_c + 0.205 ser_L_c + 0.21 pro_L_c + 0.229 asn_L_c + 0.229 asp_L_c + 0.241 thr_L_c + 0.25 gln_L_c + 0.25 glu_L_c + 0.276 ile_L_c + 0.281 arg_L_c + 0.326 lys_L_c + 0.402 val_L_c + 0.428 leu_L_c + 0.488 ala_L_c + 0.582 gly_c + 0.087 cys_L_c + 0.09 his_L_c + 0.131 tyr_L_c + 0.146 met_L_c + 0.176 phe_L_c + 0.266 g6p_c + 0.148 f6p_c + 0.052 r5p_c + 0.279 g3p_c + 0.152 3pg_c + 0.1 pep_c + 0.143 pyr_c + 4.755 accoa_c + 0.0476 akg_c + 0.0476 oaa_c + 0.121 atp_c + 0.141 gtp_c + 0.101 ctp_c + 0.102 utp_c + 0.007 ttp_c";
        ChemicalEquationSide result = ParseSubstrateEquationSide(raw_reaction, delimiters);
        REQUIRE(result.size() == 35);
        REQUIRE(result[0].name == "trp_L_c");
        REQUIRE(result[0].coefficient == Approx(0.054));
        REQUIRE(result[5].name == "thr_L_c");
        REQUIRE(result[5].coefficient == Approx(0.241));
        REQUIRE(result.back().name == "ttp_c");
        REQUIRE(result.back().coefficient == Approx(0.007));
    }

    SECTION("Without coefficients") {
        std::string raw_reaction = "a + b + 3 c + 4.0 d";
        ChemicalEquationSide  result = ParseSubstrateEquationSide(raw_reaction, delimiters);
        REQUIRE(result.size() == 4);
        REQUIRE(result[0].name == "a");
        REQUIRE(result[0].coefficient == Approx(1.0));
        REQUIRE(result[1].name == "b");
        REQUIRE(result[1].coefficient == Approx(1.0));
        REQUIRE(result[2].name == "c");
        REQUIRE(result[2].coefficient == Approx(3.0));
        REQUIRE(result[3].name == "d");
        REQUIRE(result[3].coefficient == Approx(4.0));
    }

    SECTION("Two coefficients") {
        std::string raw_reaction = "5.0 a + 228 1337 b";
        REQUIRE_THROWS_WITH(ParseSubstrateEquationSide(raw_reaction, delimiters), "There is reaction with two coefficient in a row!");
    }

    SECTION("Empty") {
        std::string raw_reaction;
        REQUIRE_THROWS_WITH(ParseSubstrateEquationSide(raw_reaction, delimiters), "");
    }

    SECTION("Two pluses") {
        std::string raw_reaction = "5.0 a + + b";
        REQUIRE_THROWS_WITH(ParseSubstrateEquationSide(raw_reaction, delimiters), "");
    }

    SECTION("Same substrate") {
        std::string raw_reaction = "5.0 a + 3.4 a";
        REQUIRE_THROWS_WITH(ParseSubstrateEquationSide(raw_reaction, delimiters), "");
    }
}

auto CompareChemicalEquationSide(const ChemicalEquationSide& lhs, const ChemicalEquationSide& rhs) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i].name != rhs[i].name) {
            return false;
        }
        bool is_equal = (lhs[i].coefficient - 0.01 < rhs[i].coefficient && lhs[i].coefficient + 0.01 > rhs[i].coefficient);
        if (!is_equal) {
            return false;
        }
        if (lhs[i].formula != rhs[i].formula) {
            return false;
        }
    }
    return true;
}

auto CompareChemicalEquations(const ChemicalEquation& lhs, const ChemicalEquation& rhs) {


    if (!CompareChemicalEquationSide(lhs.left, rhs.left)) {
        return false;
    }
    if (!CompareChemicalEquationSide(lhs.right, rhs.right)) {
        return false;
    }
    return true;
}


TEST_CASE("ParseChemicalEquation()") {
    Delimiters delimiters;
    delimiters.csv_delimiter = ',';
    delimiters.reaction_side_delimiter = "=";
    delimiters.substrate_delimiter = '+';



    SECTION("Normal") {
        std::string equation = "2 SUC = OAA + FADH2 + 0.5 NADH, 2 abcd = abcd + X + X";
        std::stringstream line(equation);
        auto result = ParseChemicalEquation(line, delimiters);

        ChemicalEquation should_be;
        ChemicalEquationSide& left = should_be.left;
        {
            Substrate suc{"SUC", "abcd", 2.0};
            left.push_back(suc);
        }

        ChemicalEquationSide& right = should_be.right;
        {
            Substrate oaa{"OAA", "abcd", 1.0};
            Substrate fadh{"FADH2", "X", 1.0};
            Substrate nadh{"NADH", "X", 0.5};
            right.push_back(oaa);
            right.push_back(fadh);
            right.push_back(nadh);
        }

        REQUIRE(CompareChemicalEquations(result, should_be));
    }

    SECTION("Missed space between coefficient and substrate") {
        std::string equation = "2SUC = OAA, 2 ab = ab";
        std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }

    SECTION("Missed space between coefficient and atoms") {
        std::string equation = "2 SUC = OAA, 2ab = ab";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }

    SECTION("Missed coefficient in atoms") {
        std::string equation = "2 SUC = OAA, ab = ab";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }

    SECTION("Missed coefficient in substrate") {
        std::string equation = "SUC = OAA, 2.0 ab = ab";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }

    SECTION("Two coefficient in a row") {
        std::string equation = "2 SUC = 0.5 0.5 OAA, 2ab = ab";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "There is reaction with two coefficient in a row!");
    }

    SECTION("Same substrate in one side") {
        std::string equation = "2 SUC + SUC = 0.5 OAA, 2ab + cd = ab";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }

    SECTION("Atoms are not consistent") {
        std::string equation = "A + B = C, ab + cd = aye";
    std::stringstream line(equation);
        REQUIRE_THROWS_WITH(ParseChemicalEquation(line, delimiters), "");
    }
}

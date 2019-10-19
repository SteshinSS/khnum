#include "catch/catch.hpp"

#include <stdio.h>
#include <iostream>
#include <sstream>

#include "parser/open_flux_parser/open_flux_utills.h"

using namespace khnum;
using namespace khnum::open_flux_parser;

const std::string folder = "../tests/open_flux_parser_test";

TEST_CASE("GetLines()", "[OpenFluxParser]") {
    const std::string subfolder = folder + "/get_lines";
    SECTION("with empty file") {
        std::vector<std::string> result = GetLines(subfolder + "/empty_file.txt");
        std::vector<std::string> empty_vector = {};
        REQUIRE(result == empty_vector);
    }
    SECTION("with two lines") {
        std::vector<std::string> result = GetLines(subfolder + "/file_with_two_lines.txt");
        std::vector<std::string> lines = {"Hello", "World!"};
        REQUIRE(result == lines);
    }
    SECTION("with next line symbol") {
        std::vector<std::string> result = GetLines(subfolder + "/file_with_newline.txt");
        std::vector<std::string> lines = {"Hello", "World!"};
        REQUIRE(result == lines);
    }
    SECTION("with empty string") {
        std::vector<std::string> result = GetLines(subfolder + "/file_with_empty_line.txt");
        std::vector<std::string> lines = {"Hello", "World!"};
        REQUIRE(result == lines);
    }
    SECTION("with empty symbols") {
        std::vector<std::string> result = GetLines(subfolder + "/file_with_empty_symbols.txt");
        std::vector<std::string> line = {"Hello   world hi"};
        REQUIRE(result == line);
    }
    SECTION("with non-existing file") {
        REQUIRE_THROWS(GetLines(subfolder + "/unexisting"));
    }
}

TEST_CASE("ParseMeasuredIsotopes()", "[OpenFluxParser]") {
    SECTION("with normal file") {
        std::vector<Emu> result = ParseMeasuredIsotopes({"Emu1:111", "OtherEmu2:1"});
        REQUIRE(result[0].name == "Emu1");
        REQUIRE(result[0].atom_states.size() == 3);
        for (size_t i = 0; i < 3; ++i) {
            REQUIRE(result[0].atom_states[i] == 1);
        }

        REQUIRE(result[1].name == "OtherEmu2");
        REQUIRE(result[1].atom_states.size() == 1);
        for (size_t i = 0; i < 3; ++i) {
            REQUIRE(result[1].atom_states[i] == 1);
        }
    }

    SECTION("with repeated emus") {
        std::vector<std::string> measured_isotopes = {"VALX:11", "ASP:1111", "VALX:1"};
        REQUIRE_THROWS(ParseMeasuredIsotopes(measured_isotopes));
    }

    SECTION("with no emus") {
        std::vector<std::string> measured_isotopes{};
        REQUIRE_THROWS(ParseMeasuredIsotopes(measured_isotopes));
    }
}

TEST_CASE("ParseOneMeasuredIsotope()", "[OpenFluxParser]") {
    SECTION("with normal emu") {
        std::string line = "VALX:1111";
        Emu emu = ParseOneMeasuredIsotope(line);
        REQUIRE(emu.name == "VALX");
        REQUIRE(emu.atom_states.size() == 4);
        for (size_t i = 0; i < 4; ++i) {
            REQUIRE(emu.atom_states[i] == 1);
        }
    }

    SECTION("with wrong atom states") {
        std::string emu_with_zero = "ASP:101";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(emu_with_zero));

        std::string emu_strange = "VALX:123";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(emu_strange));
    }

    SECTION("with no atom states") {
        std::string emu_first = "VALX:";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(emu_first));

        std::string emu_second = "QWE";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(emu_second));
    }

    SECTION("with missed next line symbol") {
        std::string long_emu = "VALX:111ASP:111";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(long_emu));
    }

    SECTION("with empty isotope") {
        std::string empty = "";
        REQUIRE_THROWS(ParseOneMeasuredIsotope(empty));
    }
}

TEST_CASE("ParseMeasurements()", "[OpenFluxParser]") {
    const std::string subfolder = folder + "/parse_measurements";

    Delimiters delimiters;
    delimiters.csv_delimiter = ',';

    SECTION("normal case") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/normal_case_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        const std::vector<std::string> raw_measurements = GetLines(subfolder + "/normal_case_measurements.txt");
        std::vector<Measurement> measurements = ParseMeasurements(raw_measurements, measured_isotopes, delimiters);

        REQUIRE(measurements.size() == 2);

        Measurement& first = measurements[0];
        REQUIRE(first.emu == measured_isotopes[0]);
        auto compare_vec = [](const std::vector<double>& lhs, const std::vector<double>& rhs) {
            REQUIRE(lhs.size() == rhs.size());
            for (size_t i = 0; i < lhs.size(); ++i) {
                REQUIRE(lhs[i] == Approx( rhs[i] ));
            }
        };

        std::vector<double> first_mid = {1.0, 2.0};
        compare_vec(first.mid, first_mid);

        std::vector<double> first_error = {0.001, 0.001};
        compare_vec(first.errors, first_error);

        Measurement& second = measurements[1];
        REQUIRE(second.emu == measured_isotopes[1]);

        std::vector<double> second_mid = {10.0, 8.0, 9.13, 32.0};
        compare_vec(second.mid, second_mid);

        std::vector<double> second_errors = {0.228, 0.228, 0.228, 0.228};
        compare_vec(second.errors, second_errors);
    }

    SECTION("too little measurements") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/little_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/little_measurements.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }

    SECTION("redundant measurements") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/redundant_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/redundant_measurements.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }

    SECTION("wrong symbols") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/wrong_symbols_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/wrong_symbols_measurements.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }

    SECTION("overflow") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/overflow_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/overflow_measurements.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }

    SECTION("integer numbers") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/integer_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/integer_measurements.txt");
        std::vector<Measurement> measurements = ParseMeasurements(raw_measurements, measured_isotopes, delimiters);

        REQUIRE(measurements[0].mid[0] == Approx(1.0));
        REQUIRE(measurements[0].mid[1] == Approx(3.0));

        REQUIRE(measurements[0].errors[0] == Approx(2.0));
        REQUIRE(measurements[0].errors[1] == Approx(4.0));
    }

    SECTION("No header") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/normal_case_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/no_header.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }

    SECTION("Empty") {
        std::vector<std::string> raw_measured_isotopes = GetLines(subfolder + "/normal_case_emus.txt");
        std::vector<Emu> measured_isotopes = ParseMeasuredIsotopes(raw_measured_isotopes);

        std::vector<std::string> raw_measurements = GetLines(subfolder + "/empty.txt");
        REQUIRE_THROWS(ParseMeasurements(raw_measurements, measured_isotopes, delimiters));
    }
}

TEST_CASE("GetCell()", "[OpenFluxParser]") {
    Delimiters delimiter;
    delimiter.csv_delimiter = ',';

    SECTION("normal line") {
        std::string string_line = "a,hello, hello world! ,   ,";
        std::stringstream line(string_line);

        std::string a = GetCell(line, delimiter);
        REQUIRE(a == "a");

        std::string hello = GetCell(line, delimiter);
        REQUIRE(hello == "hello");

        std::string hello_world = GetCell(line, delimiter);
        REQUIRE(hello_world == " hello world!");

        std::string spaces = GetCell(line, delimiter);
        REQUIRE(spaces == "   ");

        std::string empty = GetCell(line, delimiter);
        REQUIRE(empty.empty());
    }
}

TEST_CASE("ParseInputSubstrate()", "[OpenFluxParser]") {
    Delimiters delimiters;
    delimiters.csv_delimiter = ',';

    const std::string subfolder = folder + "/parse_input_substrate";
    SECTION("normal use") {
        std::vector<std::string> raw_input_substrates = GetLines(subfolder + "/normal_case.txt");
        std::vector<InputSubstrate> substrates = ParseInputSubstrates(raw_input_substrates, delimiters);

        REQUIRE(substrates.size() == 2);

        InputSubstrate& first = substrates[0];
        REQUIRE(first.name == "PYR_EX");
        REQUIRE(first.mixtures.size() == 2);
        REQUIRE(first.mixtures[0].ratio == Approx(0.5));

        auto compare_vec = [](const std::vector<double>& lhs, const std::vector<double>& rhs) {
            REQUIRE(lhs.size() == rhs.size());
            for (size_t i = 0; i < lhs.size(); ++i) {
                REQUIRE(lhs[i] == Approx(rhs[i]));
            }
        };

        std::vector<double> first_first_pat = {0.99, 0.0107, 0.0107};
        compare_vec(first.mixtures[0].fractions, first_first_pat);

        std::vector<double> first_second_pat = {0.99, 0.99, 0.99};
        compare_vec(first.mixtures[1].fractions, first_second_pat);

        InputSubstrate& second = substrates[1];
        REQUIRE(second.name == "GLU_EX");
        REQUIRE(second.mixtures.size() == 1);
        REQUIRE(second.mixtures[0].ratio == Approx(1));

        std::vector<double> second_pat = {0.99, 0.0107, 0.0107, 0.0107, 0.0107};
        compare_vec(second.mixtures[0].fractions, second_pat);
    }

    SECTION("Denormalized ratios") {
        std::vector<std::string> raw_input_substrates = GetLines(subfolder + "/denormalized.txt");
        std::vector<InputSubstrate> substrates = ParseInputSubstrates(raw_input_substrates, delimiters);

        REQUIRE(substrates[0].mixtures[0].ratio == Approx(0.5));
        REQUIRE(substrates[0].mixtures[0].ratio == Approx(0.5));
    }

    SECTION("Wrong commas") {
        std::vector<std::string> raw_input_substrates = GetLines(subfolder + "/wrong_commas.txt");
        REQUIRE_THROWS(ParseInputSubstrates(raw_input_substrates, delimiters));
    }

    SECTION("Empty file") {
        std::vector<std::string> raw_input_substrates = GetLines(subfolder + "/empty.txt");
        REQUIRE_THROWS(ParseInputSubstrates(raw_input_substrates, delimiters));
    }

    SECTION("No header") {
        std::vector<std::string> raw_input_substrates = GetLines(subfolder + "/no_header.txt");
        REQUIRE_THROWS(ParseInputSubstrates(raw_input_substrates, delimiters));
    }
}
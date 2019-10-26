#include <modeller/calculate_input_mid.h>
#include "catch/catch.hpp"
#include "modeller/modeller.h"

#include "modeller/create_emu_reactions.h"

using namespace khnum;
using namespace khnum::modelling_utills;

TEST_CASE("CalculateOneMid()", "[Modelling Utils]") {
    SECTION("One Mix, Two Masses") {
        Mixture mix;
        mix.ratio = 1.0;
        mix.fractions = {0.1, 0.3};
        InputSubstrate input_substrate;
        input_substrate.name = "my_emu";
        input_substrate.mixtures.push_back(mix);

        Emu emu;
        emu.name = "my_emu";
        emu.atom_states = {1, 1};

        auto result = CalculateOneMid(input_substrate, emu);
        Mid should_be = {0.63, 0.43, 0.03};
        REQUIRE(result.mid.size() == should_be.size());
        for (int i = 0; i < should_be.size(); ++i) {
            REQUIRE(result.mid[i] == Approx(should_be[i]));
        }
    }

    SECTION("Two Mixes, Two Masses") {
        Mixture mix;
        mix.ratio = 0.3333;
        mix.fractions = {0.1, 0.3};

        Mixture second_mix;
        second_mix.ratio = 0.6666;
        second_mix.fractions = {0.6, 0.1};

        InputSubstrate input_substrate;
        input_substrate.name = "my_emu";
        input_substrate.mixtures.push_back(mix);
        input_substrate.mixtures.push_back(second_mix);

        Emu emu;
        emu.name = "my_emu";
        emu.atom_states = {0, 1};

        auto result = CalculateOneMid(input_substrate, emu);
        Mid should_be = {0.83, 0.16};
        REQUIRE(result.mid.size() == should_be.size());
        for (int i = 0; i < should_be.size(); ++i) {
            REQUIRE(result.mid[i] == Approx(should_be[i]));
        }
    }
}

TEST_CASE("CreateOneEmuReaction", "[Modelling Utils]") {
    SECTION("Normal use") {
        Reaction reaction;
        reaction.id = 5;
        ChemicalEquation equation; // "A + 2.0 B = 1.5 C + D, ab + cd = bda + c"
        Substrate A;
        A.name = "A";
        A.formula = "ab";
        A.coefficient = 1.0;

        Substrate B;
        B.name = "B";
        B.formula = "cd";
        B.coefficient = 2.0;

        Substrate C;
        C.name = "C";
        C.formula = "bda";
        C.coefficient = 1.5;

        Substrate D;
        D.name = "D";
        D.formula = "c";

        equation.left.push_back(A);
        equation.left.push_back(B);
        equation.right.push_back(C);
        equation.right.push_back(D);

        reaction.chemical_equation = equation;

        Emu C110;
        C110.name = "C";
        C110.atom_states = {1, 0, 1};

        EmuReaction result = CreateOneEmuReaction(reaction, C, C110);
        EmuReaction should_be;
        should_be.id = reaction.id;
        should_be.right.emu = C110;

        EmuSubstrate A01;
        A01.coefficient = 1.0;
        A01.emu.name = "A";
        A01.emu.atom_states = {0, 1};
        should_be.left.push_back(A01);

        EmuSubstrate B01;
        B01.coefficient = 1.5;
        B01.emu.name = "B";
        B01.emu.atom_states = {0, 1};
        should_be.left.push_back(B01);

        REQUIRE(result == should_be);
    }
}
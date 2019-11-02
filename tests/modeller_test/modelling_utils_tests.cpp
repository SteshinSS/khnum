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
        Mid should_be = {0.63, 0.34, 0.03};
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
            REQUIRE(result.mid[i] == Approx(should_be[i]).epsilon(0.1));
        }
    }
}

TEST_CASE("CreateOneEmuReaction", "[Modelling Utils]") {
    SECTION("Complicated use") {
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
        D.coefficient = 1.0;

        equation.left.push_back(A);
        equation.left.push_back(B);
        equation.right.push_back(C);
        equation.right.push_back(D);

        reaction.chemical_equation = equation;

        Emu C110;
        C110.name = "C";
        C110.atom_states = {1, 1, 0};

        EmuReaction result = CreateOneEmuReaction(reaction, C, C110);
        EmuReaction should_be; // A:01 + 2.0 B:01 = 1.5 C:110
        should_be.id = reaction.id;
        should_be.right.emu = C110;
        should_be.right.coefficient = 1.5;

        EmuSubstrate A01;
        A01.coefficient = 1.0;
        A01.emu.name = "A";
        A01.emu.atom_states = {0, 1};
        should_be.left.push_back(A01);

        EmuSubstrate B01;
        B01.coefficient = 2.0;
        B01.emu.name = "B";
        B01.emu.atom_states = {0, 1};
        should_be.left.push_back(B01);

        REQUIRE(result == should_be);
    }

    SECTION("Convolution of three") {
        Reaction reaction;
        reaction.id = 228;
        ChemicalEquation equation; // A + B + C + D = E + F, ab + cd + ef + g = abcdef + g
        Substrate A;
        A.name = "A";
        A.formula = "ab";
        A.coefficient = 1.0;
        equation.left.push_back(A);

        Substrate B;
        B.name = "B";
        B.formula = "cd";
        B.coefficient = 1.0;
        equation.left.push_back(B);

        Substrate C;
        C.name = "C";
        C.formula = "ef";
        C.coefficient = 1.0;
        equation.left.push_back(C);

        Substrate D;
        D.name = "D";
        D.formula = "g";
        D.coefficient = 1.0;
        equation.left.push_back(D);

        Substrate E;
        E.name = "E";
        E.formula = "abcdef";
        E.coefficient = 1.0;
        equation.right.push_back(E);

        Substrate F;
        F.name = "F";
        F.formula = "g";
        F.coefficient = 1.0;
        equation.right.push_back(F);

        reaction.chemical_equation = equation;

        Emu E100101;
        E100101.name = "E";
        E100101.atom_states = {1, 0, 0, 1, 0, 1};

        EmuReaction result = CreateOneEmuReaction(reaction, E, E100101);
        EmuReaction should_be;
        should_be.id = 228;

        EmuSubstrate A10;
        A10.coefficient = 1.0;
        A10.emu.name = "A";
        A10.emu.atom_states = {1, 0};
        should_be.left.push_back(A10);

        EmuSubstrate B01;
        B01.coefficient = 1.0;
        B01.emu.name = "B";
        B01.emu.atom_states = {0, 1};
        should_be.left.push_back(B01);

        EmuSubstrate C01;
        C01.coefficient = 1.0;
        C01.emu.name = "C";
        C01.emu.atom_states = {0, 1};
        should_be.left.push_back(C01);

        should_be.right.coefficient = 1.0;
        should_be.right.emu = E100101;

        REQUIRE(result == should_be);
    }
}

TEST_CASE("CreateNewEmuReactions()", "[Modelling Utils]") {
    SECTION("Left symmetry first") {
        Reaction reaction;
        reaction.id = 1337;
        ChemicalEquation equation; // Fum + Fum = OAC, 0.5 abcd + 0.5 dcba = abcd
        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.left.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.left.push_back(Fum_second);

        Substrate OAC;
        OAC.coefficient = 1.0;
        OAC.name = "OAC";
        OAC.formula = "abcd";
        equation.right.push_back(OAC);
        reaction.chemical_equation = equation;

        Emu OAC1100;
        OAC1100.name = "OAC";
        OAC1100.atom_states = {1, 1, 0, 0};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, OAC1100);
        REQUIRE(result.size() == 2);

        EmuReaction should_be_first;
        should_be_first.id = 1337;

        EmuSubstrate Fum1100;
        Fum1100.coefficient = 0.5;
        Fum1100.emu.name = "Fum";
        Fum1100.emu.atom_states = {1, 1, 0, 0};
        should_be_first.left.push_back(Fum1100);

        should_be_first.right.emu = OAC1100;
        should_be_first.right.coefficient = 1.0;

        bool is_first_ok = (result[0] == should_be_first) || (result[1] == should_be_first);

        EmuReaction should_be_second;
        should_be_second.id = 1337;

        EmuSubstrate Fum0011;
        Fum0011.coefficient = 0.5;
        Fum0011.emu.name = "Fum";
        Fum0011.emu.atom_states = {0, 0, 1, 1};
        should_be_second.left.push_back(Fum0011);

        should_be_first.right.emu = OAC1100;
        should_be_first.right.coefficient = 1.0;

        bool is_second_ok = (result[0] == should_be_second) || (result[1] == should_be_second);
        REQUIRE( (is_first_ok && is_second_ok) );
    }

    SECTION("Left symmetry second") {
        Reaction reaction;
        reaction.id = 1337;
        ChemicalEquation equation; // Fum + Fum = OAC, 0.5 abcd + 0.5 dcba = abcd
        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.left.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.left.push_back(Fum_second);

        Substrate OAC;
        OAC.coefficient = 1.0;
        OAC.name = "OAC";
        OAC.formula = "abcd";
        equation.right.push_back(OAC);
        reaction.chemical_equation = equation;

        Emu OAC1001;
        OAC1001.name = "OAC";
        OAC1001.atom_states = {1, 0, 0, 1};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, OAC1001);
        REQUIRE(result.size() == 1);
        EmuReaction reaction_result = result.front();

        EmuReaction should_be;
        should_be.id = 1337;

        EmuSubstrate Fum1001;
        Fum1001.coefficient = 1.0;
        Fum1001.emu.name = "Fum";
        Fum1001.emu.atom_states = {1, 0, 0, 1};
        should_be.left.push_back(Fum1001);

        should_be.right.emu = OAC1001;
        should_be.right.coefficient = 1.0;

        REQUIRE(reaction_result == should_be);
    }

    SECTION("Symmetry right first") {
        Reaction reaction;
        reaction.id = 666;
        ChemicalEquation equation; // OAC = Fum + Fum, abcd = 0.5 abcd + 0.5 dcba

        Substrate OAC;
        OAC.coefficient = 1.0;
        OAC.name = "OAC";
        OAC.formula = "abcd";
        equation.left.push_back(OAC);

        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.right.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.right.push_back(Fum_second);
        reaction.chemical_equation = equation;

        Emu Fum1100;
        Fum1100.name = "Fum";
        Fum1100.atom_states = {1, 1, 0, 0};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, Fum1100);
        REQUIRE(result.size() == 2);

        EmuReaction should_be_first; // OAC1100 = Fum1100
        should_be_first.id = 666;

        EmuSubstrate OAC1100;
        OAC1100.coefficient = 1.0;
        OAC1100.emu.name = "OAC";
        OAC1100.emu.atom_states = {1, 1, 0, 0};

        should_be_first.left.push_back(OAC1100);

        should_be_first.right.coefficient = 0.5;
        should_be_first.right.emu = Fum1100;

        bool is_first_ok = (result[0] == should_be_first) || (result[1] == should_be_first);

        EmuReaction should_be_second; // OAC0011 = Fum1100
        should_be_second.id = 666;

        EmuSubstrate OAC0011;
        OAC0011.coefficient = 1.0;
        OAC0011.emu.name = "OAC";
        OAC0011.emu.atom_states = {0, 0, 1, 1};

        should_be_second.left.push_back(OAC0011);

        should_be_second.right.coefficient = 0.5;
        should_be_second.right.emu = Fum1100;

        bool is_second_ok = (result[0] == should_be_second) || (result[1] == should_be_second);

        REQUIRE((is_first_ok && is_second_ok));
    }

    SECTION("Symmetry right second") {
        Reaction reaction;
        reaction.id = 666;
        ChemicalEquation equation; // OAC = Fum + Fum, abcd = 0.5 abcd + 0.5 dcba

        Substrate OAC;
        OAC.coefficient = 1.0;
        OAC.name = "OAC";
        OAC.formula = "abcd";
        equation.left.push_back(OAC);

        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.right.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.right.push_back(Fum_second);
        reaction.chemical_equation = equation;

        Emu Fum1001;
        Fum1001.name = "Fum";
        Fum1001.atom_states = {1, 0, 0, 1};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, Fum1001);
        REQUIRE(result.size() == 1);

        EmuReaction should_be; // OAC1001 = Fum1001
        should_be.id = 666;

        EmuSubstrate OAC1001;
        OAC1001.coefficient = 1.0;
        OAC1001.emu.name = "OAC";
        OAC1001.emu.atom_states = {1, 0, 0, 1};

        should_be.left.push_back(OAC1001);

        should_be.right.coefficient = 1.0;
        should_be.right.emu = Fum1001;

        REQUIRE(result.front() == should_be);
    }

    SECTION("Symmetry first") {
        Reaction reaction;
        reaction.id = 1337;
        ChemicalEquation equation; // Suc + Suc = Fum + Fum, 0.5 abcd + 0.5 dcba = 0.5 abcd + 0.5 dcba

        Substrate Suc_first;
        Suc_first.coefficient = 0.5;
        Suc_first.name = "Suc";
        Suc_first.formula = "abcd";
        equation.left.push_back(Suc_first);

        Substrate Suc_second;
        Suc_second.coefficient = 0.5;
        Suc_second.name = "Suc";
        Suc_second.formula = "dcba";
        equation.left.push_back(Suc_second);

        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.right.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.right.push_back(Fum_second);

        Emu Fum1100;
        Fum1100.name = "Fum";
        Fum1100.atom_states = {1, 1, 0, 0};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, Fum1100);
        // 1.0 Suc1100 = 1.0 Fum1100
        // 1.0 Suc0011 = 1.0 Fum1100

        REQUIRE(result.size() == 2);

        EmuReaction should_be_first;
        should_be_first.id = 1337;

        EmuSubstrate Suc1100;
        Suc1100.coefficient = 1.0;
        Suc1100.emu.name = "Suc";
        Suc1100.emu.atom_states = {1, 1, 0, 0};
        should_be_first.left.push_back(Suc1100);

        should_be_first.right.coefficient = 1.0;
        should_be_first.right.emu = Fum1100;

        bool is_first_ok = (result[0] == should_be_first) || (result[1] == should_be_first);

        EmuReaction should_be_second;
        should_be_second.id = 1337;

        EmuSubstrate Suc0011;
        Suc0011.coefficient = 1.0;
        Suc0011.emu.name = "Suc";
        Suc0011.emu.atom_states = {0, 0, 1, 1};
        should_be_second.left.push_back(Suc0011);

        should_be_second.right.coefficient = 1.0;
        should_be_second.right.emu = Fum1100;

        bool is_second_ok = (result[0] == should_be_second) || (result[1] == should_be_second);

        REQUIRE( (is_first_ok && is_second_ok));
    }

    SECTION("Symmetry second") {
        Reaction reaction;
        reaction.id = 1337;
        ChemicalEquation equation; // Suc + Suc = Fum + Fum, 0.5 abcd + 0.5 dcba = 0.5 abcd + 0.5 dcba

        Substrate Suc_first;
        Suc_first.coefficient = 0.5;
        Suc_first.name = "Suc";
        Suc_first.formula = "abcd";
        equation.left.push_back(Suc_first);

        Substrate Suc_second;
        Suc_second.coefficient = 0.5;
        Suc_second.name = "Suc";
        Suc_second.formula = "dcba";
        equation.left.push_back(Suc_second);

        Substrate Fum_first;
        Fum_first.coefficient = 0.5;
        Fum_first.name = "Fum";
        Fum_first.formula = "abcd";
        equation.right.push_back(Fum_first);

        Substrate Fum_second;
        Fum_second.coefficient = 0.5;
        Fum_second.name = "Fum";
        Fum_second.formula = "dcba";
        equation.right.push_back(Fum_second);

        Emu Fum1001;
        Fum1001.name = "Fum";
        Fum1001.atom_states = {1, 0, 0, 1};

        std::vector<EmuReaction> result = CreateNewEmuReactions(reaction, Fum1001);

        // 2.0 Suc1001 = 1.0 Fum1001
        REQUIRE(result.size() == 1);

        EmuReaction should_be;
        should_be.id = 1337;

        EmuSubstrate Suc1001;
        Suc1001.coefficient = 2.0;
        Suc1001.emu.name = "Suc";
        Suc1001.emu.atom_states = {1, 0, 0, 1};
        should_be.left.push_back(Suc1001);

        should_be.right.coefficient = 1.0;
        should_be.right.emu = Fum1001;

        REQUIRE(result.front() == should_be);
    }
}
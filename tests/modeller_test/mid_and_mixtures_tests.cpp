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

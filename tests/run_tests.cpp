#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"
#include "utilities/debug_utills/debug_prints.h"
#include "utilities/emu.h"

TEST_CASE("First") {
    khnum::PrintEmu(khnum::Emu{"hi", {'1', '1'}});
    REQUIRE(true);
}
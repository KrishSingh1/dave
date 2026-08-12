#include "catch.hpp"
#include "guess.h"
#include <iostream>

TEST_CASE("Fluxes stores daily per-PFT fluxes independently", "[fluxes]") {
    date.init(1);
    npft = 3;
    nst = 1;

    pftlist.killall();
    for (int i = 0; i < npft; ++i) {
        Pft p = pftlist.createobj();
        p.id = i;
    }

    Soiltype soil;
    Gridcell gc;
    Stand stand(0, &gc, soil, FOREST, 1);
    Patch& patch = stand.getobj();

    date.day = 10;
    date.month = 0;
    patch.fluxes.report_flux(Fluxes::NPP, 0, 1.0);
    patch.fluxes.report_flux(Fluxes::NPP, 1, 2.0);
    patch.fluxes.report_flux(Fluxes::GPP, 1, 3.0);

    date.day = 11;
    patch.fluxes.report_flux(Fluxes::NPP, 1, 4.0);

    date.day = 364;
    patch.fluxes.report_flux(Fluxes::NPP, 0, 5.0);

    // Test boundary condition - last day, last PFT, last flux type
    Fluxes::PerPFTFluxType flux_type = Fluxes::PerPFTFluxType(Fluxes::NPERPFTFLUXTYPES - 1);
    patch.fluxes.report_flux(flux_type, npft - 1, 6.0);

    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 10, 0) == Approx(1.0));
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 10, 1) == Approx(2.0));
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::GPP, 10, 1) == Approx(3.0));
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 11, 1) == Approx(4.0));
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 10) == Approx(3.0));
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 364, 0) == Approx(5.0));
    REQUIRE(patch.fluxes.get_daily_flux(flux_type, 364, npft - 1) == Approx(6.0));

    REQUIRE(patch.fluxes.get_annual_flux(Fluxes::NPP) == Approx(12.0));
    REQUIRE(patch.fluxes.get_annual_flux(flux_type) == Approx(6.0));

    patch.fluxes.reset();
    REQUIRE(patch.fluxes.get_daily_flux(Fluxes::NPP, 10, 1) == Approx(0.0));
}
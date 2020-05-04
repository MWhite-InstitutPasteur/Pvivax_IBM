/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include <Rcpp.h>
#include "seasonality.h"

context("Seasonality") {
    test_that("the fourier switch works") {
        params p;
        p.a0 = 0.284596;
        p.a_seasonality = {
            -0.317878,
            -0.0017527,
            0.116455
        };
        p.b_seasonality = {
            -0.331361,
            0.293128,
            -0.0617547
        };
        p.dry_seas[0] = .5;
        p.kappa_seas[0] = 2.;
        p.t_peak_seas[0] = 0;
        p.use_fourier = false;
        CATCH_CHECK(seasonality(p, 0, 0) == Approx(2.666667));
        CATCH_CHECK(seasonality(p, 0, 1000) == Approx(1.62505));
        p.use_fourier = true;
        CATCH_CHECK(seasonality(p, 0, 0) == Approx(0.081423));
        CATCH_CHECK(seasonality(p, 0, 1000) == Approx(0.6370574981));
    }
}

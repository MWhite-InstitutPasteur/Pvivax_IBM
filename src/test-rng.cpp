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
#include "rng.h"
#include <Rcpp.h>

void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

context("Random number generation") {
    test_that("Exponential generation is reproducible") {
        set_seed(42);
        CATCH_REQUIRE(genexp(1) == Approx(0.1983368));
    }

    test_that("Normal generation is reproducible") {
        set_seed(42);
        CATCH_REQUIRE(gennor(0, 1) == Approx(1.370958));
    }

    test_that("Uniform generation is reproducible") {
        set_seed(42);
        CATCH_REQUIRE(genunf(0, 1) == Approx(0.914806));
    }

    test_that("Multivariate normal generation is reproducible") {
        set_seed(42);
        float covm[] = {1., 0, 0, 1.};
        float meanv[] = {0., 0.};
        long p = 2L;
        float param[10];
        setgmn(meanv, covm, p, param);
        float result[2];
        float work[2];
        genmn(param, result, work);
        CATCH_REQUIRE(result[0] == Approx(1.37096));
        CATCH_REQUIRE(result[1] == Approx(-0.564698));
    }
}

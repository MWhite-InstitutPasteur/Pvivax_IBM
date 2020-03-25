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
#include "randlib.h"
#include <Rcpp.h>

using namespace Rcpp;

void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

bool almost_equal(double a, double b, double delta) {
	return a < b + delta && a > b - delta;
}

context("Random number generation") {
  test_that("Exponential generation is reproducible") {
	  set_seed(42);
	  expect_true(almost_equal(rgenexp(1), 0.1983368, .0001));
  }

  test_that("Normal generation is reproducible") {
  	  set_seed(42);
  	  expect_true(almost_equal(rgennor(0, 1), 1.370958, .0001));
  }

  test_that("Uniform generation is reproducible") {
	  set_seed(42);
	  expect_true(almost_equal(rgenunf(0, 1), 0.914806, .0001));
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
  	  rgenmn(param, result, work);
  	  expect_true(almost_equal(result[0], 1.37096, .0001));
  	  expect_true(almost_equal(result[1], -0.564698, .0001));
  }
}

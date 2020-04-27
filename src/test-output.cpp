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
#include "output.h"

bool all_sug(Rcpp::LogicalVector x) {
    return Rcpp::is_true(Rcpp::all(x));
}

void expect_dataframes_equal(Rcpp::DataFrame observed, Rcpp::DataFrame expected) {
    auto observed_names = static_cast<Rcpp::CharacterVector>(observed.names());
    auto expected_names = static_cast<Rcpp::CharacterVector>(expected.names());
    CATCH_REQUIRE(observed_names.length() == expected_names.length());
    CATCH_REQUIRE(setequal(observed_names, expected_names));
    for (auto name : expected_names) {
        auto string_name = Rcpp::as<std::string>(name);
        Rcpp::NumericVector observed_vector = observed[string_name];
        Rcpp::NumericVector expected_vector = expected[string_name];
        CATCH_CHECK(all_sug(observed_vector == expected_vector));
    }
}


context("Model outputs") {

    test_that("Summary disaggregates population correctly") {
        simulation sim;
        population pop;
        pop.people = {
            individual(52., 1.),
            individual(29., 1.),
            individual(1352., 1.),
            individual(235., 1.),
            individual(4262., 1.)
        };
        pop.prev_groups = {
            std::pair<int,int>(-1, -1),
            std::pair<int,int>(0, 5),
            std::pair<int,int>(5, 15)
        };
        pop.incidence_groups = {
            std::pair<int,int>(-1, -1),
            std::pair<int,int>(2, 10),
        };

        pop.N_pop = 5;
        POP_summary(&pop, &sim);
        CATCH_CHECK(pop.prev_summaries[0][0] == 5);
        CATCH_CHECK(pop.prev_summaries[1][0] == 4);
        CATCH_CHECK(pop.prev_summaries[2][0] == 1);
    }

    test_that("Default model output has correct columns") {
        simulation sim;
        sim.t_vec = {365, 730};
        sim.yH_t = {
            {1, 2, 1, 1, 0, 0},
            {1, 2, 1, 1, 0, 0}
        };
        sim.yM_t = {
            {
                {1, 2, 1, 1, 0, 0},
                {1, 2, 1, 1, 0, 0},
                {1, 2, 1, 1, 0, 0}
            }, //t = 0
            {
                {1, 2, 1, 1, 0, 0},
                {1, 2, 1, 1, 0, 0},
                {1, 2, 1, 1, 0, 0}
            } //t = 1
        };
        sim.prev_groups = {
            std::pair<int,int>(-1, -1),
            std::pair<int,int>(0, 5),
            std::pair<int,int>(5, 15)
        };
        sim.incidence_groups = {
            std::pair<int,int>(-1, -1),
            std::pair<int,int>(2, 10),
        };
        sim.prev_summaries = {
            {
                {5, 2, 2, 2, 2, 2},
                {5, 4, 4, 4, 4, 4}
            },
            {
                {4, 2, 2, 2, 2, 2},
                {4, 4, 4, 4, 4, 4}
            },
            {
                {1, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0}
            }
        };
        sim.incidence_summaries = {
            {
                {0, 0, 0, 0, 0},
                {2, 2, 2, 2, 2}
            },
            {
                {0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0}
            }
        };
        sim.EIR_t = { 2.3, 5 };
        sim.LLIN_cov_t = { 4, 1 };
        sim.IRS_cov_t = { 4, 1 };
        sim.ACT_treat_t = { 4, 1 };
        sim.PQ_treat_t = { 4, 1 };
        sim.pregnant_t = { 4, 1 };
        sim.A_par_mean_t = { .3, .5 };
        sim.A_clin_mean_t = { .3, .5 };

        auto expected_frame = create_wide_dataframe({
            column_descriptor("time", Rcpp::NumericVector::create(365, 730)),
            column_descriptor("S", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("I_PCR", Rcpp::NumericVector::create(2, 2)),
            column_descriptor("I_LM", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("D", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("T", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("P", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("EL_M_far", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("LL_M_far", Rcpp::NumericVector::create(2, 2)),
            column_descriptor("P_M_far", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("S_M_far", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("E_M_far", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("I_M_far", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("EL_M_pun", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("LL_M_pun", Rcpp::NumericVector::create(2, 2)),
            column_descriptor("P_M_pun", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("S_M_pun", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("E_M_pun", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("I_M_pun", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("EL_M_kol", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("LL_M_kol", Rcpp::NumericVector::create(2, 2)),
            column_descriptor("P_M_kol", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("S_M_kol", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("E_M_kol", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("I_M_kol", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("N_pop", Rcpp::NumericVector::create(5, 5)),
            column_descriptor("PvPR_PCR", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvPR_LM", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("Pv_clin", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvHR", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvHR_batch", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("new_PCR", Rcpp::NumericVector::create(0, 2)),
            column_descriptor("new_LM", Rcpp::NumericVector::create(0, 2)),
            column_descriptor("new_D", Rcpp::NumericVector::create(0, 2)),
            column_descriptor("new_T", Rcpp::NumericVector::create(0, 2)),
            column_descriptor("N_pop_0_5", Rcpp::NumericVector::create(4, 4)),
            column_descriptor("PvPR_PCR_0_5", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvPR_LM_0_5", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("Pv_clin_0_5", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvHR_0_5", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("PvHR_batch_0_5", Rcpp::NumericVector::create(2, 4)),
            column_descriptor("N_pop_5_15", Rcpp::NumericVector::create(1, 1)),
            column_descriptor("PvPR_PCR_5_15", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("PvPR_LM_5_15", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("Pv_clin_5_15", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("PvHR_5_15", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("PvHR_batch_5_15", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("new_PCR_2_10", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("new_LM_2_10", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("new_D_2_10", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("new_T_2_10", Rcpp::NumericVector::create(0, 0)),
            column_descriptor("EIR", Rcpp::NumericVector::create(2.3, 5)),
            column_descriptor("LLIN_cov", Rcpp::NumericVector::create(4, 1)),
            column_descriptor("IRS_cov", Rcpp::NumericVector::create(4, 1)),
            column_descriptor("ACT_treat", Rcpp::NumericVector::create(4, 1)),
            column_descriptor("PQ_treat", Rcpp::NumericVector::create(4, 1)),
            column_descriptor("pregnant", Rcpp::NumericVector::create(4, 1)),
            column_descriptor("A_par", Rcpp::NumericVector::create(.3, .5)),
            column_descriptor("A_clin", Rcpp::NumericVector::create(.3, .5))
        });

        expect_dataframes_equal(create_output_frame(sim), expected_frame);
    }
}

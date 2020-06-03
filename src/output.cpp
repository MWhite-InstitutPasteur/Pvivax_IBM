/*
 * output.cpp
 *
 *  Created on: 6 Apr 2020
 *      Author: gc1610
 */

#include "output.h"

Rcpp::DataFrame create_output_frame(const simulation& sim) {
    auto columns = std::vector<column_descriptor>();
    columns.push_back(column_descriptor("time", Rcpp::wrap(sim.t_vec)));
    const auto human_states = std::vector<std::string>{
        "S",
        "I_PCR",
        "I_LM",
        "D",
        "T",
        "P"
    };

    for (auto i = 0u; i < human_states.size(); ++i) {
        auto column = std::vector<double>();
        for (auto t = 0u; t < sim.yH_t.size(); ++t) {
            column.push_back(sim.yH_t[t][i]);
        }
        columns.push_back(column_descriptor(human_states[i], Rcpp::wrap(column)));
    }

    const auto mosquito_varieties = std::vector<std::string>{
        "far",
        "pun",
        "kol"
    };
    const auto mosquito_states = std::vector<std::string>{
        "EL",
        "LL",
        "P",
        "S",
        "E",
        "I"
    };

    for (auto i = 0u; i < mosquito_varieties.size(); ++i) {
        for (auto prev_type = 0u; prev_type < mosquito_states.size(); ++prev_type) {
            auto column = std::vector<double>();
            std::stringstream colname;
            colname << mosquito_states[prev_type] << "_M_"  << mosquito_varieties[i];
            for (auto t = 0u; t < sim.yM_t.size(); ++t) {
                column.push_back(sim.yM_t[t][i][prev_type]);
            }
            columns.push_back(
                column_descriptor(
                    colname.str(),
                    Rcpp::wrap(column)
                )
            );
        }
    }

    const auto prev_types = std::vector<std::vector<std::string>>{
        {
            "N_pop",
            "PvPR_PCR",
            "PvPR_LM",
            "Pv_clin",
            "PvHR",
            "PvHR_batch"
        },
        {
            "new_PCR",
            "new_LM",
            "new_D",
            "new_T"
        }
    };

    const auto prev_groups = std::vector<std::vector<std::pair<int,int>>>{
        sim.prev_groups,
        sim.incidence_groups
    };

    const auto prev_summaries = std::vector<std::vector<std::vector<std::vector<int>>>>{
        sim.prev_summaries,
        sim.incidence_summaries
    };

    for (auto inc_prev = 0u; inc_prev < prev_groups.size(); ++inc_prev) {
        for (auto prev_group = 0u; prev_group < prev_groups[inc_prev].size(); ++prev_group) {
            for (auto prev_type = 0u; prev_type < prev_types[inc_prev].size(); ++prev_type) {
                auto column = std::vector<double>();
                std::stringstream colname;
                if (prev_group == 0) {
                    colname << prev_types[inc_prev][prev_type];
                } else {
                    colname << prev_types[inc_prev][prev_type];
                    colname << '_' << prev_groups[inc_prev][prev_group].first;
                    colname << '_' << prev_groups[inc_prev][prev_group].second;
                }
                const auto & summary = prev_summaries[inc_prev];
                for (auto t = 0u; t < summary.size(); ++t) {
                    column.push_back(summary[t][prev_group][prev_type]);
                }
                columns.push_back(
                    column_descriptor(
                        colname.str(),
                        Rcpp::wrap(column)
                    )
                );
            }
        }
    }

    columns.push_back(column_descriptor("EIR", Rcpp::wrap(sim.EIR_t)));
    columns.push_back(column_descriptor("LLIN_cov", Rcpp::wrap(sim.LLIN_cov_t)));
    columns.push_back(column_descriptor("IRS_cov", Rcpp::wrap(sim.IRS_cov_t)));
    columns.push_back(column_descriptor("ACT_treat", Rcpp::wrap(sim.ACT_treat_t)));
    columns.push_back(column_descriptor("PQ_treat", Rcpp::wrap(sim.PQ_treat_t)));
    columns.push_back(column_descriptor("pregnant", Rcpp::wrap(sim.pregnant_t)));
    columns.push_back(column_descriptor("A_par", Rcpp::wrap(sim.A_par_mean_t)));
    columns.push_back(column_descriptor("A_clin", Rcpp::wrap(sim.A_clin_mean_t)));

    return create_wide_dataframe(columns);
}

Rcpp::DataFrame create_wide_dataframe(std::vector<column_descriptor> columns) {
    auto df = Rcpp::DataFrame();
    for (const auto& column : columns) {
        df.push_back(column.second, column.first);
    }
    return df;
}

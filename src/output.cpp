/*
 * output.cpp
 *
 *  Created on: 6 Apr 2020
 *      Author: gc1610
 */

#include "output.h"

Rcpp::DataFrame create_output_frame(simulation* sim) {
    auto columns = std::vector<column_descriptor>();
    columns.push_back(column_descriptor("time", Rcpp::wrap(sim->t_vec)));
    auto human_states = vector<std::string>{
        "S",
        "I_PCR",
        "I_LM",
        "D",
        "T",
        "P"
    };
    auto human_state_columns = std::vector<std::vector<double>>(human_states.size());
    for (auto t = 0; t < sim->yH_t.size(); ++t) {
        for (auto i = 0; i < human_states.size(); ++i) {
            human_state_columns[i].push_back(sim->yH_t[t][i]);
        }
    }
    for (auto i = 0; i < human_states.size(); ++i) {
        columns.push_back(column_descriptor(human_states[i], Rcpp::wrap(human_state_columns[i])));
    }
    return create_wide_dataframe(columns);
}

Rcpp::DataFrame create_wide_dataframe(std::vector<column_descriptor> columns) {
    auto df = Rcpp::DataFrame();
    for (const auto& column : columns) {
        df.push_back(column.second, column.first);
    }
    return df;
}

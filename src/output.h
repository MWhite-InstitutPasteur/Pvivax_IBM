/*
 * output.h
 *
 *  Created on: 8 Apr 2020
 *      Author: gc1610
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_

#include "model.h"

using column_descriptor = std::pair<std::string, Rcpp::NumericVector>;
using suffix_vector_pair = std::pair<std::string, std::vector<std::vector<int>>>;

Rcpp::DataFrame create_output_frame(const simulation&);
Rcpp::DataFrame create_wide_dataframe(std::vector<column_descriptor>);


#endif /* SRC_OUTPUT_H_ */

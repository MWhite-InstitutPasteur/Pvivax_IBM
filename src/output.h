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

Rcpp::DataFrame create_output_frame(simulation*);
Rcpp::DataFrame create_wide_dataframe(std::vector<column_descriptor>);


#endif /* SRC_OUTPUT_H_ */

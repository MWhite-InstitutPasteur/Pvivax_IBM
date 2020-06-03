/*
 * seasonality.cpp
 *
 *  Created on: 30 Apr 2020
 *      Author: gc1610
 */

#include "model.h"
#include <math.h>
#include <Rcpp.h>


double seasonality(params& p, size_t g, size_t t) {
    if (p.use_fourier) {
        double theta = 0.01721421 * (t - p.t_peak_seas[g]);
        double y = p.a0;
        for(auto i = 0u; i < p.a_seasonality.size(); ++i) {
            y += p.a_seasonality[i] * cos(theta*(i+1)) +
                p.b_seasonality[i] * sin(theta*(i+1));
        }
        return y;
    }

    auto denominator = exp(
        gammln(0.5) + gammln(p.kappa_seas[g] + 0.5
    ) - gammln(p.kappa_seas[g] + 1.0)) / 3.14159265359;

    auto numerator = p.dry_seas[g] +
        (1 - p.dry_seas[g]) *
        pow(
            0.5 + 0.5 * cos(0.01721421 * (t - p.t_peak_seas[g])),
            p.kappa_seas[g]
        );

    return numerator / denominator;
}

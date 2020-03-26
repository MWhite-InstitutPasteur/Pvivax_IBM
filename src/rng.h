/*
 * rng.h
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#ifndef SRC_RNG_H_
#define SRC_RNG_H_

double genexp(double av);
double gennor(double av, double sd);
double genunf(double low, double high);
void genmn(float *parm,float *x,float *work);
void setgmn(float *meanv,float *covm,long p,float *parm);

#endif /* SRC_RNG_H_ */

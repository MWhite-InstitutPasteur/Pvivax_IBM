/*
 * rng.h
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#ifndef SRC_RNG_H_
#define SRC_RNG_H_

float genexp(float av);
float gennor(float av,float sd);
float genunf(float low,float high);
void genmn(float *parm,float *x,float *work);
void setgmn(float *meanv,float *covm,long p,float *parm);

#endif /* SRC_RNG_H_ */

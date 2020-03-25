/*
 * rng.h
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#ifndef SRC_RNG_H_
#define SRC_RNG_H_

float rgenexp(float av);
float rgennor(float av,float sd);
float rgenunf(float low,float high);
void rgenmn(float *parm,float *x,float *work);

#endif /* SRC_RNG_H_ */

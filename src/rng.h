/*
 * rng.h
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#ifndef SRC_RNG_H_
#define SRC_RNG_H_

extern float genexp(float av);
extern void genmn(float *parm,float *x,float *work);
extern float gennor(float av,float sd);
extern float genunf(float low,float high);
extern void setall(long iseed1,long iseed2);
extern void setgmn(float *meanv,float *covm,long p,float *parm);

#endif /* SRC_RNG_H_ */

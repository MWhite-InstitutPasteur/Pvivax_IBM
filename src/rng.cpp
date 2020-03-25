/*
 * rng.cpp
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#include <Rcpp.h>

float rgenexp(float av) {
	return Rcpp::rexp(1, av)[0];
}

float rgennor(float av, float sd) {
	return Rcpp::rnorm(1, av, sd)[0];
}

float rgenunf(float low, float high) {
	return Rcpp::runif(1, low, high)[0];
}


// Adapted from randlib.c
void rgenmn(float *parm,float *x,float *work) {
static long i,icount,j,p,D1,D2,D3,D4;
static float ae;

    p = (long) (*parm);
/*
     Generate P independent normal deviates - WORK ~ N(0,1)
*/
    for(i=1; i<=p; i++) *(work+i-1) = rgennor(0, 1);
    for(i=1,D3=1,D4=(p-i+D3)/D3; D4>0; D4--,i+=D3) {
/*
     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
      decomposition of the desired covariance matrix.
          trans(A)(1,1) = PARM(P+2)
          trans(A)(2,1) = PARM(P+3)
          trans(A)(2,2) = PARM(P+2+P)
          trans(A)(3,1) = PARM(P+4)
          trans(A)(3,2) = PARM(P+3+P)
          trans(A)(3,3) = PARM(P+2-1+2P)  ...
     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
*/
        icount = 0;
        ae = 0.0;
        for(j=1,D1=1,D2=(i-j+D1)/D1; D2>0; D2--,j+=D1) {
            icount += (j-1);
            ae += (*(parm+i+(j-1)*p-icount+p)**(work+j-1));
        }
        *(x+i-1) = ae+*(parm+i);
    }
}

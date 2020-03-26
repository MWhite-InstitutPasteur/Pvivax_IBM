/*
 * rng.cpp
 *
 *  Created on: 25 Mar 2020
 *      Author: gc1610
 */

#include <Rcpp.h>

float genexp(float av) {
	return R::rexp(av);
}

float gennor(float av, float sd) {
	return R::rnorm(av, sd);
}

float genunf(float low, float high) {
	return R::runif(low, high);
}


void setgmn(float *meanv,float *covm,long p,float *parm)
/*
**********************************************************************
     void setgmn(float *meanv,float *covm,long p,float *parm)
            SET Generate Multivariate Normal random deviate
                              Function
      Places P, MEANV, and the Cholesky factoriztion of COVM
      in GENMN.
                              Arguments
     meanv --> Mean vector of multivariate normal distribution.
     covm   <--> (Input) Covariance   matrix    of  the  multivariate
                 normal distribution
                 (Output) Destroyed on output
     p     --> Dimension of the normal, or length of MEANV.
     parm <-- Array of parameters needed to generate multivariate norma
                deviates (P, MEANV and Cholesky decomposition of
                COVM).
                1 : 1                - P
                2 : P + 1            - MEANV
                P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
               Needed dimension is (p*(p+3)/2 + 1)
**********************************************************************
*/
{
extern void spofa(float *a,long lda,long n,long *info);
static long T1;
static long i,icount,info,j,D2,D3,D4,D5;
    T1 = p*(p+3)/2+1;
/*
     TEST THE INPUT
*/
    if(!(p <= 0)) goto S10;
    Rcpp::Rcerr << "Value of P: " << p << std::endl;
    Rcpp::stop("P nonpositive in SETGMN");
S10:
    *parm = p;
/*
     PUT P AND MEANV INTO PARM
*/
    for(i=2,D2=1,D3=(p+1-i+D2)/D2; D3>0; D3--,i+=D2) *(parm+i-1) = *(meanv+i-2);
/*
      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
*/
    spofa(covm,p,p,&info);
    if(!(info != 0)) goto S30;
    Rcpp::stop("COVM not positive definite in SETGMN");
S30:
    icount = p+1;
/*
     PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
          COVM(1,1) = PARM(P+2)
          COVM(1,2) = PARM(P+3)
                    :
          COVM(1,P) = PARM(2P+1)
          COVM(2,2) = PARM(2P+2)  ...
*/
    for(i=1,D4=1,D5=(p-i+D4)/D4; D5>0; D5--,i+=D4) {
        for(j=i-1; j<p; j++) {
            icount += 1;
            *(parm+icount-1) = *(covm+i-1+j*p);
        }
    }
}

// Adapted from randlib.c
void genmn(float *parm,float *x,float *work) {
/*
**********************************************************************
	 void genmn(float *parm,float *x,float *work)
			  GENerate Multivariate Normal random deviate
							  Arguments
	 parm --> Parameters needed to generate multivariate normal
			   deviates (MEANV and Cholesky decomposition of
			   COVM). Set by a previous call to SETGMN.
			   1 : 1                - size of deviate, P
			   2 : P + 1            - mean vector
			   P+2 : P*(P+3)/2 + 1  - upper half of cholesky
									   decomposition of cov matrix
	 x    <-- Vector deviate generated.
	 work <--> Scratch array
							  Method
	 1) Generate P independent standard normal deviates - Ei ~ N(0,1)
	 2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
	 3) trans(A)E + MEANV ~ N(MEANV,COVM)
**********************************************************************
*/
static long i,icount,j,p,D1,D2,D3,D4;
static float ae;

    p = (long) (*parm);
/*
     Generate P independent normal deviates - WORK ~ N(0,1)
*/
    for(i=1; i<=p; i++) *(work+i-1) = gennor(0, 1);
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

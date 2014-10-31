/*********************************************************************** 
 * preconditioning related functions.
 * Note: there exist already polynomial preconditioning operators
 *       related to the chebychev polynomials in poly_precon.h files
 *       that was developed earlier by Carsten
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/

#ifndef _PRECON_H
#define _PRECON_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "linalg_eo.h"
#include "memalloc.h"
#include "start.h"

void cheb_poly_precon(spinor * const R, spinor * const S, matrix_mult f, const int N, const double evmin, const double evmax, const int k);
/*
 R = T_{k}(Q)*S (k>=0) where S is input spinor and Q is the operator given by the matrix-vector multipliction where Q*v is given by the
 matrix-vector multiplication opertor f. T_{k}(Q) is the chebychev polynomial.
  
 For -1 =< t =< 1, the Chebyshev polynomials are given by

       T_0(t) = 1    
       T_1(t) = t
       T_j(t) = 2*t*T_{j-1}(t) - T_{j-2}(t) where j=2,3,...

and the roots of T_n(t) in the interval [-1,1] for n>=1 are given by

       t_i = cos(pi/2 *(2*i-1)/j) for i=1,2,..,j   and j>=1

In our case, we choose an interval [evmin,evmax] which will require shifting the polynomial.

So we define t = 2/(evmax-evmin)*x - (evmax+evmin)/(evmax-evmin)

where x here is our matrix Q.

So, as a polynomial of Q, we have 
    T_0(Q) = 1

    T_1(Q) = 2/(evmax-evmin)*Q - (evmax+evmin)/(evmax-evmin)

    T_j(Q) = 4/(evmax-evmin)*Q*T_{j-1}(Q) - 2*(evmax+evmin)/(evmax-evmin)*T_{j-1}(Q) - T_{j-2}(Q)  


and the roots of this polynomial are 

    q_i = (evmax+evmin)/2+ (evmax-evmin)/2*cos(pi/2 *(2*i-1)/j)   for i=1,2,..,j and j>=1 


R: output spinor
S: input spinor
f: matrix-vector multiplication operator
N: size of the spinor
evmin: minimum value of the interval
evmax: maximum value of the interval
k: order of the Chebeychev polynomial
*/

#endif

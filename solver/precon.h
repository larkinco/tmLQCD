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
 R = T_{k+1}(Q)*S (k>=0) where S is input spinor and Q is the operator given by the matrix-vector multipliction where Q*v is given by the
 matrix-vector multiplication opertor av. T_{k+1}(Q) is the chebychev polynomial given by:
  
       evmin and evmax are minimum and maximum eigenvalues of the operator under consideration (assumed to be Hermitian here)  

       theta = (evmax+evmin)/2,   delta=(evmax-evmin)/2
       sigma_0 =1 ,     sigma_1=theta/delta
       rho_{-1}=0 ,     rho_0=1/(2*sigma_1)
       rho_k   = 1/(2*sigma_1 - rho_{k-1})

       T_{-1}(t) = 0,    T_0(t)=1,  T_1(t)=1-t/theta
       
       T_{k+1}(t) = rho_k*{ 2*(sigma_1-t/delta)*T_k(t) - rho_{k-1}*T_{k-1}) 
*/


void cheb_poly_precon_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double evmin, const double evmax, const int k);
//R = T_{k+1}(Q)*(Q*S). Thi is the preconditioned Dirac operator.


#endif

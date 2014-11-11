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

void cheb_poly_precon_residual(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k);
/*
This is the polynomial T_{k}(x) = 1 -xP_{k-1}(x) defined such that |1-xP_{k-1}(x)| is minimized over the interval [a,b]. P_{k-1}(x) is the approximate
inverse over this interval and T(x) is the residual polynomial. When solving a linear system Ay=b with initial residual r_0=b-Ay_0, we see that the 
resiudal of the preconditioned system P(A)A y = P(A)b has a residual T(A)r_0. For eigenvalue problems, we can use T_{k}(A) as a filter to minimize
the contribution coming from eigenvalues in the interval [a,b] and enhance the other part of the spectrum. This will be used with ARPACK for example.
For the linear system on the other hand we will use P(A) as a preconditioner and it is defined in the function cheb_poly_precon_op. Note that the roots
of T(x) are related to the roots of the Chebyshev polynomial of the first kind.

Propertires of the Chebyshev polynomials of the first kind C_k(x) for x in the interval [-1,1] and k=0,1,2,3,...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

k=0: C_0(x) =1
k=1: C_1(x) = x
...
three-term recurrence:   C_{k+1}(x) = 2xC_k(x) - C_{k-1}

C_k(x) has k simple roots in the interval [-1,1] given by
x_l = cos(pi/2 (2l-1)/k) for k=1,2,3,... and l=1,2,..,k

In the interval [a,b] we can define a shifted polynomial by C_k(d(x) where d(x) = 2x/(b-a) - (b+a)/(b-a)
such that d(a)=-1, and d(b)=+1. 

d(x) = x/delta - theta/delta  with delta = (b-a)/2 and theta=(b+a)/2

R_k(x) = C_k(d(x) is the shifted Chebyshev polynomial with 

R_0(x) =1
R_1(x) = d(x)

R_{k+1}(x) = 2d(x)R_k(x)-R_{k-1}(x)

T_k(x) is the normalized shifted Chebyshev polynomial given by:

For k=1,2,3,....   T_k(x) = R_k(x)/sigma_k(x)

with sigma_0 = 1, sigma_1=d(0) = - theta/delta
sigma_{k+1} =  -2theta/delta sigma_k - sigma_{k-1}

These relations define the polynomial T_k(x)

The roots of T_k(x) are the same as the roots of C_k(d(x)) and are given by

x_l = (b-a)/2*[cos(pi/2 (2*l-1)/k)+(b+a)/(b-a)]

*/


//This is the version of the preconditioner to be used with linear system solution
void cheb_poly_precon_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k);
/*
Computes the action of the polynomial P_k(x) which minimizes the |1-xP_k(x)| over the interval [a,b]. This polynomial is used as 
an approximation to the inverse of the matrix over that interval and can be used as a preconditioner for CG. The recurrence relations
could be derived similar to the residual polynomial above and are given by:

sigma_0 = 1, sigma_1=d(0) = - theta/delta
sigma_{k+1} =  -2theta/delta sigma_k - sigma_{k-1}

P_0(Q) = -1/theta

P_1(Q) = - (4*theta+2*Q)/(2*theta^2-delta^2)

P_{k+1}(Q) = -(2*sigma_{k+1}/(delta*sigma_{k+2}))- 2*sigma_{k+1}/(delta*sigma_{k+2})*(theta-Q)*P_k(Q) - sigma_k/sigma_{k+2}*P_{k-1}(Q)

*/




void cheb_poly_roots(_Complex double *roots, const int k, const double a, const double b);
/*
roots of the shifted Chebyshev polynomial in the interval [a,b]
*/


#endif

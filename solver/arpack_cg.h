/*****************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Deflating CG using eigenvectors computed using ARPACK
 *
 * Author: A.M. Abdel-Rehim (amabdelrehim@gmail.com)
 *         Novemebr 2014
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************************/
/* A sample input is given in the sample-input folder */

#ifndef _ARPACK_CG_H
#define _ARPACK_CG_H

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/eigenvalues_arpack.h"


int arpack_cg(

     //solver params
     const int N,             /*(IN) Number of lattice sites for this process*/
     const int nrhs,          /*(IN) Number of right-hand sides to be solved*/ 
     const int nrhs1,         /*(IN) First number of right-hand sides to be solved using tolerance eps_sq1*/ 
     spinor * const x,        /*(IN/OUT) initial guess on input, solution on output for this RHS*/
     spinor * const b,        /*(IN) right-hand side*/
     matrix_mult f,           /*(IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply*/
     const double eps_sq1,    /*(IN) squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
     const double eps_sq,     /*(IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
     const double res_eps_sq, /*(IN) suqared tolerance for restarting cg */
     const int rel_prec,      /*(IN) 0 for using absoute error for convergence
                                     1 for using relative error for convergence*/
     const int maxit,         /*(IN) Maximum allowed number of iterations to solution for the linear system*/

     //parameters for arpack
     const int nev,                 /*(IN) number of eigenvectors to be computed by arpack*/
     const int ncv,                 /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
     double arpack_eig_tol,         /*(IN) tolerance for computing eigenvalues with arpack */
     int arpack_eig_maxiter,        /*(IN) maximum number of iterations to be used by arpack*/
     int kind,                      /*(IN) 0 for eigenvalues with smallest real part "SR"
                                           1 for eigenvalues with largest  real part "LR"
                                           2 for eigenvalues with smallest absolute value "SM"
                                           3 for eigenvalues with largest absolute value  "LM"
                                           4 for eigenvalues with smallest imaginary part "SI"
                                           5 for eigenvalues with largest imaginary part  "LI"*/
     int comp_evecs,                /*(IN) 0 no computation of the residuals of eigenvectors of the operator
                                           1 compute residuals of the eigenevtors of the operator*/
     int acc,                       /*(IN) 0 no polynomial acceleration
                                           1 use polynomial acceleration*/
     int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
     double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     int read_basis,               /*(IN) 0 compute deflation basis using arpack, 1 read deflation basis from disk */
     int store_basis,              /*(IN) when using arpack to compute eigenvectors use this  
                                          option to store basis vectors to disk such that they can be read later
                                          0 don't store basis vectors
                                          1 store basis vectors */
     char *basis_fname,            /*(IN)file name used to read/store the basis vectors
                                         file names will be of the format
                                         basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
     int basis_prec,               /*(IN)precision used to write the basis vectors
                                         0 single precision
                                         1 double precision*/
     char *arpack_logfile,           /*(IN) file name to be used for printing out debugging information from arpack*/
     int top_suscept_exact_flag
     );


//solver to compute eigenvectors only
int arpack(
     //solver params
     const int N,             /*(IN) Number of lattice sites for this process*/
     matrix_mult f,           /*(IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply*/
     //parameters for arpack
     const int nev,                 /*(IN) number of eigenvectors to be computed by arpack*/
     const int ncv,                 /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
     double arpack_eig_tol,         /*(IN) tolerance for computing eigenvalues with arpack */
     int arpack_eig_maxiter,        /*(IN) maximum number of iterations to be used by arpack*/
     int kind,                      /*(IN) 0 for eigenvalues with smallest real part "SR"
                                           1 for eigenvalues with largest  real part "LR"
                                           2 for eigenvalues with smallest absolute value "SM"
                                           3 for eigenvalues with largest absolute value  "LM"
                                           4 for eigenvalues with smallest imaginary part "SI"
                                           5 for eigenvalues with largest imaginary part  "LI"*/
     int acc,                       /*(IN) 0 no polynomial acceleration
                                           1 use polynomial acceleration*/
     int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
     double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     int store_basis,              /*(IN) when using arpack to compute eigenvectors use this  
                                          option to store basis vectors to disk such that they can be read later
                                          0 don't store basis vectors
                                          1 store basis vectors */
     char *basis_fname,            /*(IN)file name used to read/store the basis vectors
                                         file names will be of the format
                                         basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
     int basis_prec,               /*(IN)precision used to write the basis vectors
                                         0 single precision
                                         1 double precision*/
     char *arpack_logfile,           /*(IN) file name to be used for printing out debugging information from arpack*/
     int top_suscept_exact_flag,
	char *v_g5_logfile

	);

#endif 

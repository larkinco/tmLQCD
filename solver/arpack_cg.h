/*****************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * Deflating CG using eigenvectors computed using ARPACK
 * eigenvectors used correspond to those with smallest magnitude
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
     const int ncv,                 /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv =< 2*nev */
     double arpack_eig_tol,         /*(IN) tolerance for computing eigenvalues with arpack */
     int arpack_eig_maxiter,        /*(IN) maximum number of iterations to be used by arpack*/
     int kind,                      /*(IN) 0 for smallest eigenvalues, 1 for largest eigenvalues*/
     int acc,                       /*(IN) option for using polynomial acceleration:
                                           0 no acceleration
                                           1 compute the eigenvectors of the acceleration polynomial T_k(f)
                                           2 use roots of chebyshev polynomial as shifts but compute eigenvalues of f*/
     int cheb_k,                    /*(IN) degree of the Chebyshev polynomial when acc=1 (irrelevant if acc=0 and arpack_initresid=0)*/
     int emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0&arpack_initresid=0)
                                           lower bound such that eigenvalues in the interval [emin,emax] will be damped
                                           while eigenvalues outside this interval will be enhanced */
     int emax,                       /*(IN) upper end of the interval where polynomial acceleration will be used.
                                           Follow the same cases as emin above.*/
     int arpack_initresid           /*(IN) 0 means arpack will start with some random vector when computing eigenvectors
                                           1 means a starting vector will be provided by the user. In this case a vector
                                             built of of the T_k(f)v will be used where v is a random vector and T_k is
                                             the chebyshev polynomial used for acceleration.*/ 
     );

#endif 

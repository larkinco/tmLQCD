/*********************************************************************** 
 * Interface for solving the eigenvalue problem A*x=lambda*x 
 * for a complex matrix A using ARPACK and PARPACK. The matrix
 * is accessed through a matrix-vector multiplication function.
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * For reference see the driver programs zndrv1 in the EXAMPLES 
 * subdriectories of ARPACK and PARPACK.
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/

#ifndef _EIGENVALUES_ARPACK
#define _EIGENVALUES_ARPACK

#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "linalg/arpack.h"
#include "memalloc.h"
#include "solver/precon.h"


void evals_arpack(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, matrix_mult av, int *info, int *nconv);
/*
  compute nev eigenvectors using ARPACK and PARPACK
  n     : size of the lattice
  nev   : number of eigenvectors requested.
  ncv   : size of the subspace used to compute eigenvectors (nev+1) =< ncv 
  which : which eigenvectors to compute. Choices are:
          LM: largest magnitude
          SM: smallest magnitude
          LA: largest real component
          SA: smallest real compoent
          LI: largest imaginary component
          SI: smallest imaginary component
  evals : Computed eigenvalues. Size is ncv complex doubles.
  v     : Computed eigenvectors. Size is ncv*n spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 



//this version computes the eigenvectors of a shifted operator A+lambda*I
void evals_arpack_shift(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, matrix_mult av, int *info, int *nconv);
/*
  compute nev eigenvectors using ARPACK and PARPACK
  n     : size of the lattice
  nev   : number of eigenvectors requested.
  ncv   : size of the subspace used to compute eigenvectors (nev+1) =< ncv 
  which : which eigenvectors to compute. Choices are:
          LM: largest magnitude
          SM: smallest magnitude
          LA: largest real component
          SA: smallest real compoent
          LI: largest imaginary component
          SI: smallest imaginary component
  evals : Computed eigenvalues. Size is ncv complex doubles.
  v     : Computed eigenvectors. Size is ncv*n spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 





void evals_arpack_poly_hermitian(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, 
                                 matrix_mult av, double evmin, double evmax, int cheb_k, int *info, int *nconv);
/*
  compute nev eigenvectors for the Hermitian operator with polynomial preconditioning using Ceyshev polynomials using ARPACK and PARPACK 
  n     : size of the lattice
  nev   : number of eigenvectors requested.
  ncv   : size of the subspace used to compute eigenvectors (nev+1) =< ncv 
  which : which eigenvectors to compute. Choices are:
          LM: largest magnitude
          SM: smallest magnitude
          LA: largest real component
          SA: smallest real compoent
          LI: largest imaginary component
          SI: smallest imaginary component
  evals : Computed eigenvalues. Size is ncv complex doubles.
  v     : Computed eigenvectors. Size is ncv*n spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  cheb_k : the degree of the chebychev polynomial.
  evmin  : smallest eigenvalue of the chebyshev polynomial
  evmax  : largest eigenvalue of the chebyshev polynomial
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 


//In this version, polynomial acceleration is used but implemented through passing the roots of chebyshev
//polynomials to be used by arpack instead of the explicit shifts
void evals_arpack_poly_hermitian_shifts(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, 
                                 matrix_mult av, double evmin, double evmax, int *info, int *nconv);
/*
  compute nev eigenvectors for the Hermitian operator with polynomial preconditioning using Ceyshev polynomials using ARPACK and PARPACK 
  n     : size of the lattice
  nev   : number of eigenvectors requested.
  ncv   : size of the subspace used to compute eigenvectors (nev+1) =< ncv
  which : which eigenvectors to compute. Choices are:
          LM: largest magnitude
          SM: smallest magnitude
          LA: largest real component
          SA: smallest real compoent
          LI: largest imaginary component
          SI: smallest imaginary component
  evals : Computed eigenvalues. Size is ncv complex doubles.
  v     : Computed eigenvectors. Size is ncv*n spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  evmin  : smallest eigenvalue of the chebyshev polynomial
  evmax  : largest eigenvalue of the chebyshev polynomial
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 







#endif

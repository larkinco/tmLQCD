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


/*compute nev eigenvectors using ARPACK and PARPACK*/
void evals_arpack(
  int n, 
  int nev, 
  int ncv, 
  int which,
  int use_acc,
  int init_resid_arpack, 
  int cheb_k,
  double amin,
  double amax,
  _Complex double *evals, 
  spinor *v,
  _Complex double *workl,
  double tol, 
  int maxiter, 
  matrix_mult av, 
  int *info, 
  int *nconv);
/*
  compute nev eigenvectors using ARPACK and PARPACK
  n     : (IN) size of the local lattice
  nev   : (IN) number of eigenvectors requested.
  ncv   : (IN) size of the subspace used to compute eigenvectors (nev+1) =< ncv < 12*n
          where 12n is the size of the matrix under consideration
  which : (IN) which eigenvectors to compute. Choices are:
          0: smallest magnitude
          1: largest magintude
  use_acc: (IN) specify the polynomial acceleration mode
                0 no acceleration
                1 use acceleration by computing the eigenvectors of a shifted-normalized chebyshev polynomial
                2 use acceleration by using the roots of Chebyshev polynomial as shifts
  init_resid_arpack: (IN) specify the initial residual passed to arpack for computing eigenvectors
                          0 arpack uses a random intiial vector
                          1 provide a starting vector for arpack which is obtained by chebyshev
                            polynomial in order to enhance the components of the requested eiegenvectors
  cheb_k: (IN) degree of the chebyshev polynomial to be used for acceleration (irrelevant when use_acc=0 and init_resid_arpack=0)
  amin,amax: (IN) bounds of the interval [amin,amax] for the acceleration polynomial (irrelevant when use_acc=0 and init_resid_arpack=0)
  evals : (OUT) array of size nev+1 which has the computed Ritz values
  workl : (OUT) work array that has needed information about the schur decomposition and Ritz values
                size 3*ncv^2 + 5*ncv
  v     : orthonormal basis (schur vectors) of the eigenvectors. Size is ncv*ldv (ldv includes the communication buffer) spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 

#endif

/*********************************************************************** 
 * Interface for solving the eigenvalue problem A*x=lambda*x 
 * for a complex matrix A using ARPACK and PARPACK. The matrix
 * is accessed through a matrix-vector multiplication function.
 *
 * author: A.M. Abdel-Rehim
 *
 * For reference see the driver programs zndrv1 and in the EXAMPLES 
 * subdriectories of ARPACK and PARPACK.
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/

#ifndef _EIGENVALUES_ARPACK
#define _EIGENVALUES_ARPACK
#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "linalg/fortran.h"


void * alloc_aligned_mem(size_t size);
//allocate a size memory aligned on a number 
//of bytes=ALIGN_BASE boundary



void evals_arpack(int is_eo, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, matrix_mult av);
/*
  compute nev eigenvectors using ARPACK and PARPACK
  iseo  : 1 means we are solving the even-odd preconditioned system
          0 means we are solving the full size problem.
  nev   : number of eigenvectors requested.
  ncv   : size of the subspace used to compute eigenvectors (nev+1) =< ncv =< 2*nev
  which : which eigenvectors to compute. Choices are:
          LM: largest magnitude
          SM: smallest magnitude
          LA: largest real component
          SA: smallest real compoent
          LI: largest imaginary component
          SI: smallest imaginary component
  evals : Computed eigenvalues. Size is ncv complex doubles (allocated by evals_arpack).
  v     : Computed eigenvectors. Size is ncv*ldv spinors (allocated by evals_arpack).
          If using half_spinor then ldv=VOLUME (is_eo=0) or VOLUME/2 (is_eo=1).
          If not using half_spinor then ldv=VOLUMEPLUSRAND if is_eo=0
          or (VOLUMEPLUSRAND/2 if is_eo=1).

  evals and v can be freed afterwards by the calling program if needed
 
  tol   : Requested tolerance for the accuracy of the computed eigenvectors.
          A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
*/ 


#endif

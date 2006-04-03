#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<math.h>
#include"complex.h"
#include"solver/lu_solve.h"

/* Solve M a = b by LU decomposition with partial pivoting */

double norm2(complex* b, const int Nvec) {
  int i;
  double res=0.;
  for(i = 0; i < Nvec; i++) {
    res+= b[i].re*b[i].re + b[i].im*b[i].im;
  }
  return(res);
}

void LUSolve( const int Nvec, complex * M, const int ldM, complex * b) {
  int i, j, k, maxrow, row;
  complex * b_local, * y;
  double maxnorm;
  complex tmp, sum_LU;

  b_local = (complex*)malloc(Nvec*sizeof(complex));
  y = (complex*)malloc(Nvec*sizeof(complex));
  for(i = 0; i < Nvec; i++) {
    b_local[i] = b[i];
  }

  /* -----------------------------------------------------------------
   * LU Decompose M_local, in place (Crone's algorithm?)
   * It's in Numerical Recipes but also a more understandable
   * description can be found at:
   *          http://csep10.phys.utk.edu/guidry/
   *               phys594/lectures/linear_algebra/lanotes/node3.html
   * 
   * OR look in your favourite Matrix Analysis text
   * -----------------------------------------------------------------
   *
   * -------------------------------------------------------------
   * Start LU Decomp. Definition. 0-th row of U is 0-th row of M
   *   and L_{i,i} = 1 for all i
   * 
   * So we start with the 1-th (2nd) row
   * ------------------------------------------------------------ */
    
  for(i = 1; i < Nvec; i++) { 
    
    /* ------------------------------------------------------------
     * Parital Pivot: Find the row with the largest element in the
     * ith-column and make that the i-th row. This swaps rows.
     * so I don't need to reorder the unknowns, but I do need 
     * to reorder the b_local
     * ------------------------------------------------------------*/
    maxnorm = norm2(&M[i*ldM+i], 1);
    maxrow = i;
    
    /*  Compare norms with other elements in column j for row i+1.N */
    for(row = i+1; row < Nvec; row++) {
      if ( norm2(&M[row*ldM + i],1) > maxnorm ) {
	/*  Norm of M(j,i) is bigger, store it as the maximum */
	/*  and store its index */
	maxnorm = norm2(&M[row*ldM + i],1);
	maxrow = row;
      }
    }
    
    /*  If the element with maximum norm is not in row i, swap */
    /*  its row with row i */
    if( maxrow != i ) {
      
      /*  Swap rows i and maxindex */
      for(j = 0; j < Nvec; j++ ) {
	tmp = M[i*ldM + j];
	M[i*ldM + j] = M[maxrow*ldM + j];
	M[maxrow*ldM + j] = tmp;
      }
      
      /*  Swap elems of b */
      tmp = b_local[i];
      b_local[i] = b_local[maxrow];
      b_local[maxrow] = tmp;
    }
    
    /* -------------------------------------------------------- 
     * End of pivoting code
     * -------------------------------------------------------- 
    
    
     * --------------------------------------------------------
     * Work out elements of L & U in place in M for row i
     * -------------------------------------------------------- */
    for(j = 0; j < i; j++) { 
      
      _complex_set(sum_LU, 0., 0.);
      for(k = 0; k < j; k++) {
	/*  sum_LU += M(i,k)*M(k,j); */
	_mult_assign_complex(tmp, M[i*ldM + k], M[k*ldM + j]); 
	_add_complex(sum_LU, tmp);
	
      }
      /* M(i,j) -= sum_LU; */
      _diff_complex(M[i*ldM + j], sum_LU);
      /* M(i,j) /= M(j,j); */
      _div_complex(tmp, M[i*ldM + j], M[j*ldM + j]);
      M[i*ldM + j] = tmp;
    }
    
    for(j=i; j < Nvec; j++) { 
      _complex_set(sum_LU, 0., 0.);
      for(k = 0; k < i; k++) {
	_mult_assign_complex(tmp, M[i*ldM + k], M[k*ldM + j]); 
	_add_complex(sum_LU, tmp); 
      }
      /* M(i,j) -= sum_LU; */
      _diff_complex(M[i*ldM+j], sum_LU);
    }
  }
  
  /* ----------------------------------------------------
   * LU Decomp finished. M now holds the 
   *   U matrix in its diagonal and superdiagonal elements
   *   and the subdiagonal elements of the L matrix in its
   *   subdiagonal. Recall that the Diagonal elements of L 
   *   are chosen to be 1
   * -----------------------------------------------------
   
   * Solve L y = b by forward substitution */
  y[0] = b[0];
  for(i = 1; i < Nvec; i++) { 
    y[i] = b_local[i];
    for(j = 0; j < i; j++) { 
      _diff_assign_complex(y[i], M[i*ldM+j],y[j]);
    }
  }
  
  /*  Solve U a = y by back substitution */
  /* a[Nvec-1] = y[Nvec-1] / M(Nvec-1, Nvec-1); */
  _div_complex(tmp, y[Nvec-1], M[(Nvec-1)*ldM + (Nvec-1)]);
  b[Nvec-1] = tmp;
  
  for(i = Nvec-2; i >= 0; i--) { 
    tmp = y[i];
    for(j = i+1; j < Nvec; j++) { 
      /* tmp -= M(i,j)*b[j]; */
      _diff_assign_complex(tmp, M[i*ldM + j], b[j]);
    }
    _div_complex(b[i], tmp, M[i*ldM+i]);
  }
  free(b_local);
  free(y);
}


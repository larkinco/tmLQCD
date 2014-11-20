/*****************************************************************************
 * Copyright (C) 2014 Abdou M. Abdel-Rehim
 *
 * Deflating CG using eigenvectors computed using ARPACK
 * eigenvectors used correspond to those with smallest magnitude
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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef MPI
# include <mpi.h>
#endif

#include "global.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver_field.h"
#include "solver/arpack_cg.h"



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
     double arpack_eig_tol,                /*(IN) tolerance for computing eigenvalues with arpack */
     int arpack_eig_maxiter,                /*(IN) maximum number of iterations to be used by arpack*/
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
     )
{ 

  //Static variables and arrays.
  static int ncurRHS=0;                  /* current number of the system being solved */                   
  static spinor *evecs,*ax,*r;                  
  static _Complex double *evals,*H,*HU,*initwork;
  int *IPIV; 
  static int info_arpack=0;
  static int nconv_arpack=0; //number of converged eigenvectors as returned by arpack

  int i,j;

  int parallel;        /* for parallel processing of the scalar products */
  #ifdef MPI
    parallel=1;
  #else
    parallel=0;
  #endif


  _Complex double c1,c2;

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2; 

  double et1,et2;  //timers for computing eigenvetors using arpack

  //before solving 
  if(ncurRHS==0){ 
    //call arpack
    et1=gettime();
    evals = (_Complex double *) alloc_aligned_mem(ncv*sizeof(_Complex double));
    evecs = (spinor *) alloc_aligned_mem(ncv*N*sizeof(spinor));
    evals_arpack(N,nev,ncv,kind,acc,arpack_initresid,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv_arpack);
    et2=gettime();

    if(info_arpack != 0){ //arpack didn't converge
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"WARNING: ARPACK didn't converge. No deflation will be done\n");
    }
    
    if(info_arpack == 0){
      if(g_proc_id == g_stdio_proc)
        fprintf(stdout,"WARNING: ARPACK has computed %d eigenvectors\n",nconv_arpack);
      H        = (_Complex double*) alloc_aligned_mem(nconv_arpack*nconv_arpack*sizeof(_Complex double ));
      HU       = (_Complex double*) alloc_aligned_mem(nconv_arpack*nconv_arpack*sizeof(_Complex double ));
      IPIV     = (int *) alloc_aligned_mem(nconv_arpack*sizeof(int));
      initwork = (_Complex double *) alloc_aligned_mem(nconv_arpack*sizeof(_Complex double));
      if(g_proc_id == g_stdio_proc)
        fprintf(stdout,"ARPACK time: %-f\n",et2-et1);
    }
    ax     = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
    r      = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));

    //compute the elements of the hermitian matrices H and HU
    for(i=0; i<nconv_arpack; i++)
    {

       assign(r,&evecs[i*N],N);
       f(ax,r);
       for(j=i; j<nconv_arpack; j++)
       {
          c1 = scalar_prod(&evecs[j*N],ax,N,parallel);
          H[j+nconv_arpack*i] = c1;
          H[i+nconv_arpack*j] = conj(c1);
          HU[j+nconv_arpack*i] = c1;
          HU[i+nconv_arpack*j] = conj(c1);
       }
     }
  }
    
  double eps_sq_used,restart_eps_sq_used;  //tolerance squared for the linear system

  double cur_res; //current residual squared

  /*increment the RHS counter*/
  ncurRHS = ncurRHS +1; 

  //set the tolerance to be used for this right-hand side 
  if(ncurRHS > nrhs1){
    eps_sq_used = eps_sq;
  }
  else{
    eps_sq_used = eps_sq1;
  }
  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
    fprintf(stdout, "System %d, eps_sq %e\n",ncurRHS,eps_sq_used); 
    fflush(stdout);
  } 
  
  /*---------------------------------------------------------------*/
  /* Call init-CG until this right-hand side converges             */
  /*---------------------------------------------------------------*/
  double wt1,wt2,wE,wI;
  double normsq,tol_sq;
  int flag,maxit_remain,numIts,its;

  wE = 0.0; wI = 0.0;     /* Start accumulator timers */
  flag = -1;    	  /* System has not converged yet */
  maxit_remain = maxit;   /* Initialize Max and current # of iters   */
  numIts = 0;  
  restart_eps_sq_used=res_eps_sq;

  while( flag == -1 )
  {
    
    if(info_arpack==0)
    {
      /* --------------------------------------------------------- */
      /* Perform init-CG with evecs vectors                        */
      /* xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		   */
      /* --------------------------------------------------------- */
      wt1 = gettime();

      /*r0=b-Ax0*/
      f(ax,x); /*ax = A*x */
      diff(r,b,ax,N);  /* r=b-A*x */
	
      /* deflate:
         Importnat to note that here we assume that the vectors returned 
         by arpack are eigenvectors (not just an orthogonal basis).
         Also we are assuming that H=V^\dagger A V is diagonal
         For a general situation H is not diagonal and a solution of 
         a linear system is needed as in eigCG. For simplicity and
         efficiency, it is assumed that H is diagonal and that the
         arpack interface routine returns eigenvectors. This will save
         us the need to solve a small linear system for every right-hand side
         and no need to call lapack routines.
         In case we need to, we should build H after calling arpack,
         then make a cholesky factorization once. Then a back subistution
         is used for every right hand side.
      */
     
      int ONE=1;
      int info_lapack;

      /* x = x + evecs*inv(H)*evecs'*r */
      for(int i=0; i < nconv_arpack; i++)
      {
         initwork[i]= scalar_prod(&evecs[i*N],r,N,parallel);
      }

      /* solve the linear system H y = c */
      _FT(zgesv) (&nconv_arpack,&ONE,HU,&nconv_arpack,IPIV,initwork,&nconv_arpack,&info_lapack);

      if(info_lapack != 0)
      {
         if(g_proc_id == g_stdio_proc) {
            fprintf(stderr, "Error in ZGESV:, info =  %d\n",info_lapack); 
            fflush(stderr);
         }
         exit(1);
      }

      /* x = x + evecs*inv(H)*evecs'*r */
      for(i=0; i<nconv_arpack; i++)
      {
        assign_add_mul(x,&evecs[i*N],initwork[i],N);
      }

      /* compute elapsed time and add to accumulator */

      wt2 = gettime();
      wI = wI + wt2-wt1;
      
    }/* if(info_arpack == 0) */


    //which tolerance to use
    if(eps_sq_used > restart_eps_sq_used)
    {
       tol_sq = eps_sq_used;
       flag   = 1; //shouldn't restart again
    }
    else
    {
       tol_sq = restart_eps_sq_used;
    }

    wt1 = gettime();
    its = cg_her(x,b,maxit_remain,tol_sq,rel_prec,N,f); 
          
    wt2 = gettime();

    wE = wE + wt2-wt1;

    //check convergence
    if(its == -1)
    {
       //cg didn't converge
       if(g_proc_id == g_stdio_proc) {
         fprintf(stderr, "CG didn't converge within the maximum number of iterations in arpack_cg. Exiting...\n"); 
         fflush(stderr);
         exit(1);
       }
    } 
    else
    {
       numIts += its;   
       maxit_remain = maxit - numIts; //remaining number of iterations
       restart_eps_sq_used = restart_eps_sq_used*res_eps_sq; //prepare for the next restart
    }
    
  }
  /* end while (flag ==-1)               */
  
  /* ---------- */
  /* Reporting  */
  /* ---------- */
  /* compute the exact residual */
  f(ax,x); /* solver_field[0]= A*x */
  diff(r,b,ax,N);  /* solver_filed[1]=b-A*x */	
  normsq=square_norm(r,N,parallel);
  if(g_debug_level > 0 && g_proc_id == g_stdio_proc)
  {
    fprintf(stdout, "For this rhs:\n");
    fprintf(stdout, "Total initCG Wallclock : %-f\n", wI);
    fprintf(stdout, "Total cg Wallclock : %-f\n", wE);
    fprintf(stdout, "Iterations: %-d\n", numIts); 
    fprintf(stdout, "Actual Resid of LinSys  : %e\n",normsq);
  }


  //free memory if this was your last system to solve
  if(ncurRHS == nrhs){
    free(evecs);
    free(evals);
    free(ax);
    free(r);
    free(initwork);
  }


  return numIts;
}
 


      

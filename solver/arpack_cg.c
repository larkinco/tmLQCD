/*****************************************************************************
 *
 * Deflating CG using eigenvectors computed using ARPACK
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
 * Author: Abdou Abdel-Rehim
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
# include "io/spinor.h"
# include "read_input.h"

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
                                           1 for eigenvalues with largest real part "LR"
                                           2 for eigenvalues with smallest absolute value "SM"
                                           3 for eigenvalues with largest absolute value "LM"
                                           4 for eigenvalues with smallest imaginary part "SI"
                                           5 for eigenvalues with largest imaginary part  "LI"*/
     int comp_evecs,                /*(IN) 0 don't compute the eiegnvalues and their residuals of the original system 
                                           1 compute the eigenvalues and the residuals for the original system (the orthonormal baiss
                                             still be used in deflation and they are not overwritten).*/
     int acc,                       /*(IN) 0 no polynomial acceleration
                                           1 use polynomial acceleration*/
     int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
     double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     int read_basis,               /*(IN) 0 compute deflation basis using arpack, 1 read deflation basis from disk */
     int store_basis,              /*(IN) option to store basis vectors to disk such that they can be read later
                                          0 don't store basis vectors
                                          1 store basis vectors */
     char *basis_fname,            /*(IN)file name used to read/store the basis vectors
                                         file names will be of the format
                                         basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
     int basis_prec,               /*(IN)precision used to write the basis vectors
                                         0 single precision
                                         1 double precision*/
     char *arpack_logfile           /*(IN) file name to be used for printing out debugging information from arpack*/
     )
{ 

  //Static variables and arrays.
  static int  ncurRHS=0; /* current number of the system being solved */                   
  static void   *_ax,*_r,*_tmps1,*_tmps2,*_zero_spinor;                  
  static spinor *ax,*r,*tmps1,*tmps2,*zero_spinor;                  
  static _Complex double *evecs,*evals,*initwork; 
  static double *evalsA;  
  static int info_arpack=0;
  static int nconv=0; //number of converged eigenvectors as returned by arpack
  int i,j;
  complex double c1,c2;
  double d1,d2,d3;
  double et1,et2;  //timing variables
  int prec,status,rstat;
  int parallel;    /* for parallel processing of the scalar products */
  #ifdef MPI
    parallel=1;
  #else
    parallel=0;
  #endif

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2; 

  //----------------------------------------------------------------
  //if this is the first right hand side, allocate the memory needed
  //---------------------------------------------------------------- 
  if(ncurRHS==0){ 
    #if (defined SSE || defined SSE2 || defined SSE3)
    _ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_ax==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _ax inside arpack_cg.\n");
       exit(1);
    }
    else
       {ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

    _r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_r==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _r inside arpack_cg.\n");
       exit(1);
    }
    else
       {r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps1==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _tmps1 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps2==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _tmps2 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}

    _zero_spinor = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_zero_spinor==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _zero_spinor inside arpack_cg.\n");
       exit(1);
    }
    else
       {zero_spinor  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


    #else
    ax = (spinor *) malloc(LDN*sizeof(spinor));
    r  = (spinor *) malloc(LDN*sizeof(spinor));
    tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
    tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
    zero_spinor = (spinor *) malloc(LDN*sizeof(spinor));
    
    if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) || (zero_spinor==NULL) )
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for ax,r,tmps1,tmps2,zero_spinor inside arpack_cg.\n");
       exit(1);
    }
    #endif
    zero_spinor_field(zero_spinor,LDN); //this will be the even part of the eiegnvector (currently is zero)


    evecs     = (_Complex double *) malloc(ncv*12*N*sizeof(_Complex double)); //note: no extra buffer 
    evals     = (_Complex double *) malloc(ncv*sizeof(_Complex double)); 
    initwork  = (_Complex double *) malloc(nev*sizeof(_Complex double)); 
    evalsA    = (double *) malloc(nev*sizeof(double)); 

    if((evecs == NULL)  || (evals==NULL) || (evalsA==NULL))
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for evecs, evals, or evalsA inside arpack_cg.\n");
       exit(1);
    }

    //read or compute the eigenvectors
    if(read_basis)
    {
       et1=gettime();
       nconv=nev;
       for(j=0; j<nev; j++)
       {
            char filename[500],*header_type=NULL; 
            READER *reader=NULL;
            uint64_t bytes;
	    sprintf(filename, "%s.%05d", basis_fname, j);
            if( (rstat = read_spinor(zero_spinor, r, filename, 0)) != 0) {
              fprintf(stderr, "read_spinor failed with return value %d", rstat);
              exit(-7);
            }
            assign_spinor_to_complex(&evecs[j*12*N],r,N); 
       } //for(j=0;...)
       et2=gettime();
       if(g_proc_id == g_stdio_proc)
          fprintf(stdout,"Finished reading deflation basis in %e seconds\n",et2-et1); 
    }
    else //eigenvectors will be computed
    {
       et1=gettime();
       evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv,arpack_logfile);
       et2=gettime();

       if(info_arpack != 0){ //arpack didn't converge
         if(g_proc_id == g_stdio_proc)
           fprintf(stderr,"WARNING: ARPACK didn't converge. exiting..\n");
         return -1;
       }
    
       if(g_proc_id == g_stdio_proc)
       {
          fprintf(stdout,"ARPACK has computed %d eigenvectors\n",nconv);
          fprintf(stdout,"ARPACK time: %+e\n",et2-et1);
       }
       //------------------------------------------------
       //compute the eigenvalues of A and their residuals
       //------------------------------------------------
       if(g_proc_id == g_stdio_proc)
       {fprintf(stdout,"Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
        fprintf(stdout,"=============================================================\n");
        fflush(stdout);}

       for(i=0; i<nconv; i++)
       {
          assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
          f(ax,r);
          c1 = scalar_prod(r,ax,N,parallel);
          d1 = square_norm(r,N,parallel);
          evalsA[i] = creal(c1)/d1;
          mul_r(tmps1,evalsA[i],r,N);
          diff(tmps2,ax,tmps1,N);
          d2= square_norm(tmps2,N,parallel);
          d3= sqrt(d2/d1);	    
          if(g_proc_id == g_stdio_proc)
          {fprintf(stdout,"Eval[%06d]: %22.15E rnorm: %22.15E\n", i, evalsA[i], d3); fflush(stdout);}
       }

       //write the eigenvectors to disk if needed
       if(store_basis){
         for(i=0; i<nconv; i++)
         {
            assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
            WRITER *writer = NULL;
	    int append = 0;
	    char fname[256];	
	    paramsPropagatorFormat *format = NULL;
	    int precision;
	    int numb_flavs = 1;
            if(basis_prec==0)
              precision = 32;
            else
              precision = 64;

	    sprintf(fname, "%s.%05d", basis_fname, i);
	    construct_writer(&writer, fname, append);
	    format = construct_paramsPropagatorFormat(precision, numb_flavs);
	    write_propagator_format(writer, format);
	    free(format);	    
	    int status = write_spinor(writer, &zero_spinor, &r, numb_flavs, precision);
	    destruct_writer(writer);
            if(g_proc_id == g_stdio_proc)
               {fprintf(stdout,"finished writing eigenvector %d to file %s\n", i,fname); fflush(stdout);}
         } 
       } //if(store_basis)
     } //else for if(read_basis) 
  } //if(ncurRHS==0)
    
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
  int info_lapack;

  wE = 0.0; wI = 0.0;     /* Start accumulator timers */
  flag = -1;    	  /* System has not converged yet */
  maxit_remain = maxit;   /* Initialize Max and current # of iters   */
  numIts = 0;  
  restart_eps_sq_used=res_eps_sq;

  while( flag == -1 )
  {
    
    if(nconv > 0)
    {
      /* --------------------------------------------------------- */
      /* Perform init-CG with evecs vectors                        */
      /* xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		   */
      /* --------------------------------------------------------- */
      wt1 = gettime();

      /*r0=b-Ax0*/
      f(ax,x); /*ax = A*x */
      diff(r,b,ax,N);  /* r=b-A*x */

      /* x = x + evecs*inv(H)*evecs'*r */
      for(int i=0; i < nconv; i++)
      {
         assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
         initwork[i]= scalar_prod(tmps1,r,N,parallel);
         initwork[i] /= evalsA[i]; //assuming evecs are eigenvectors
         assign_add_mul(x,tmps1,initwork[i],N);
      }

      /* compute elapsed time and add to accumulator */
      wt2 = gettime();
      wI = wI + wt2-wt1;  
    }/* if(nconv > 0) */


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
  f(ax,x); /* ax= A*x */
  diff(r,b,ax,N);  /* r=b-A*x */	
  normsq=square_norm(r,N,parallel);
  if(g_debug_level > 0 && g_proc_id == g_stdio_proc)
  {
    fprintf(stdout, "For this rhs:\n");
    fprintf(stdout, "Total initCG Wallclock : %+e\n", wI);
    fprintf(stdout, "Total cg Wallclock : %+e\n", wE);
    fprintf(stdout, "Iterations: %-d\n", numIts); 
    fprintf(stdout, "Actual Resid of LinSys inside arpack_cg  : %+e\n",normsq);
  }


  //free memory if this was your last system to solve
  if(ncurRHS == nrhs){
    #if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
    free(_ax);  free(_r);  free(_tmps1); free(_tmps2);
    #else
    free(ax); free(r); free(tmps1); free(tmps2);
    #endif
    free(evecs); free(evals); free(evalsA);
    free(initwork); 
  }


  return numIts;
}
 

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
                                           1 for eigenvalues with largest real part "LR"
                                           2 for eigenvalues with smallest absolute value "SM"
                                           3 for eigenvalues with largest absolute value "LM"
                                           4 for eigenvalues with smallest imaginary part "SI"
                                           5 for eigenvalues with largest imaginary part  "LI"*/
     int acc,                       /*(IN) 0 no polynomial acceleration
                                           1 use polynomial acceleration*/
     int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
     double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
     int store_basis,              /*(IN) option to store basis vectors to disk such that they can be read later
                                          0 don't store basis vectors
                                          1 store basis vectors */
     char *basis_fname,            /*(IN)file name used to read/store the basis vectors
                                         file names will be of the format
                                         basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
     int basis_prec,               /*(IN)precision used to write the basis vectors
                                         0 single precision
                                         1 double precision*/
     char *arpack_logfile           /*(IN) file name to be used for printing out debugging information from arpack*/
     )
{


  //Static variables and arrays.
  static int  ncurRHS=0; /* current number of the system being solved */                   
  static void   *_ax,*_r,*_tmps1,*_tmps2,*_zero_spinor;                  
  static spinor *ax,*r,*tmps1,*tmps2,*zero_spinor;                  
  static _Complex double *evecs,*evals; 
  static double *evalsA;  
  static int info_arpack=0;
  static int nconv=0; //number of converged eigenvectors as returned by arpack
  int i,j;
  complex double c1,c2,c3;
  double d1,d2,d3;
  double et1,et2;  //timing variables
  int prec,status,rstat;
  int parallel;    /* for parallel processing of the scalar products */
  #ifdef MPI
    parallel=1;
  #else
    parallel=0;
  #endif

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2; 

  //----------------------------------------------------------------
  //if this is the first right hand side, allocate the memory needed
  //---------------------------------------------------------------- 
  if(ncurRHS==0){ 
    #if (defined SSE || defined SSE2 || defined SSE3)
    _ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_ax==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _ax inside arpack_cg.\n");
       exit(1);
    }
    else
       {ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

    _r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_r==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _r inside arpack_cg.\n");
       exit(1);
    }
    else
       {r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps1==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _tmps1 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps2==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _tmps2 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}

    _zero_spinor = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_zero_spinor==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for _zero_spinor inside arpack_cg.\n");
       exit(1);
    }
    else
       {zero_spinor  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


    #else
    ax = (spinor *) malloc(LDN*sizeof(spinor));
    r  = (spinor *) malloc(LDN*sizeof(spinor));
    tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
    tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
    zero_spinor = (spinor *) malloc(LDN*sizeof(spinor));
    
    if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) || (zero_spinor==NULL) )
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for ax,r,tmps1,tmps2,zero_spinor inside arpack_cg.\n");
       exit(1);
    }
    #endif
    zero_spinor_field(zero_spinor,LDN); //this will be the even part of the eiegnvector (currently is zero)


    evecs     = (_Complex double *) malloc(ncv*12*N*sizeof(_Complex double)); //note: no extra buffer 
    evals     = (_Complex double *) malloc(ncv*sizeof(_Complex double)); 
    evalsA    = (double *) malloc(nev*sizeof(double)); 

    if((evecs == NULL)  || (evals==NULL) || (evalsA==NULL))
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"insufficient memory for evecs, evals, or evalsA inside arpack_cg.\n");
       exit(1);
    }

    et1=gettime();
    evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv,arpack_logfile);
    et2=gettime();

    if(info_arpack != 0){ //arpack didn't converge
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"WARNING: ARPACK didn't converge. exiting..\n");
      return -1;
    }
    
    if(g_proc_id == g_stdio_proc)
    {
       fprintf(stdout,"ARPACK has computed %d eigenvectors\n",nconv);
       fprintf(stdout,"ARPACK time: %+e\n",et2-et1);
    }

    //------------------------------------------------
    //compute the eigenvalues of A and their residuals
    //------------------------------------------------
    if(g_proc_id == g_stdio_proc)
    {fprintf(stdout,"Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
     fprintf(stdout,"=============================================================\n");
     fflush(stdout);}

    for(i=0; i<nconv; i++)
    {
       assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
       f(ax,r);
       c1 = scalar_prod(r,ax,N,parallel);
       d1 = square_norm(r,N,parallel);
       evalsA[i] = creal(c1)/d1;
       mul_r(tmps1,evalsA[i],r,N);
       diff(tmps2,ax,tmps1,N);
       d2= square_norm(tmps2,N,parallel);
       d3= sqrt(d2/d1);	    
       if(g_proc_id == g_stdio_proc)
       {fprintf(stdout,"Eval[%06d]: %22.15E rnorm: %22.15E\n", i, evalsA[i], d3); fflush(stdout);}
    }

    //write the eigenvectors to disk if needed
    if(store_basis){
      for(i=0; i<nconv; i++)
      {
         assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
         WRITER *writer = NULL;
	 int append = 0;
	 char fname[256];	
	 paramsPropagatorFormat *format = NULL;
	 int precision;
	 int numb_flavs = 1;
         if(basis_prec==0)
           precision = 32;
         else
           precision = 64;

	 sprintf(fname, "%s.%05d", basis_fname, i);
	 construct_writer(&writer, fname, append);
	 format = construct_paramsPropagatorFormat(precision, numb_flavs);
	 write_propagator_format(writer, format);
	 free(format);	    
	 int status = write_spinor(writer, &zero_spinor, &r, numb_flavs, precision);
	 destruct_writer(writer);
         if(g_proc_id == g_stdio_proc)
            {fprintf(stdout,"finished writing eigenvector %d to file %s\n", i,fname); fflush(stdout);}
      } 
    } //if(store_basis)
    #if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
    free(_ax);  free(_r);  free(_tmps1); free(_tmps2);
    #else
    free(ax); free(r); free(tmps1); free(tmps2);
    #endif
    free(evecs); free(evals); free(evalsA);
  } //if(ncurRHS==0)
    

  /*increment the RHS counter*/
  ncurRHS = ncurRHS +1; 

  return 0;
 
}
 


      

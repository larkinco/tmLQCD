/*********************************************************************** 
 * Interface for solving the eigenvalue problem A*x=lambda*x 
 * for a complex matrix A using ARPACK and PARPACK. The matrix
 * is accessed through a matrix-vector multiplication function.
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * For reference see the driver programs zndrv1 and in the EXAMPLES 
 * subdriectories of ARPACK and PARPACK.
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/
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
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver/eigenvalues_arpack.h"


//ARPACK driver routines for computing eigenvectors 
static void _FT(znaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                          _Complex double *resid, int *ncv, _Complex double *v, int *ldv, 
                          int *iparam, int *ipntr, _Complex double *workd, _Complex double *workl, 
                          int *lworkl, double *rwork, int *info );

static void _FT(zneupd) (int *comp_evecs,char *cA, int *select, _Complex double *evals, 
                         _Complex double *v, int *ldv, _Complex double *sigma, _Complex double *workev, 
                         char *bmat, int *n, char *which, int *nev, double *tol, _Complex double *resid, 
                         int *ncv, _Complex double *v1, int *ldv1, int *iparam, int *ipntr, 
                         _Complex double *workd, _Complex double *workl, int *lworkl, double *rwork, int *ierr);

//PARPACK routines
static void _FT(pznaupd) (int *comm, int *ido, char *bmat, int *n, char *which, int *nev, double *tol, 
                         _Complex double *resid, int *ncv, _Complex double *v, int *ldv, int *iparam, 
                         int *ipntr, _Complex double *workd, _Complex double *workl, int *lworkl,
                         double *rwork, int *info );

static void _FT(pzneupd) (int *comm, int *comp_evecs,char *cA, int *select, _Complex double *evals, 
                          _Complex double *v, int *ldv, _Complex double *sigma, _Complex double *workev, 
                          char *bmat, int *n, char *which, int *nev, double *tol, _Complex double *resid, 
                          int *ncv, _Complex double *v1, int *ldv1, int *iparam, int *ipntr, 
                          _Complex double *workd, _Complex double *workl, int *lworkl, double *rwork, int *ierr);



void evals_arpack(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, matrix_mult av, int *info, int *nconv)
{
   int N,ldv;
  
   if(n==VOLUME) //full 
     ldv= VOLUMEPLUSRAND;
   else          //even-odd
     ldv= VOLUMEPLUSRAND/2;
   
   N=12*n;  //size of the matrix as complex double

   //check input
   if(nev>=N) nev=N-1;
   if(ncv < (nev+1)) ncv = nev+1;

 
   //parameters for the ARPACK routines--see documentation or source files
   int ido,lworkl,j,ierr,ishfts,mode;
   _Complex double sigma;
   int comp_evecs=1; //option to compute the eigenvectors (always set to true here)
   lworkl  = 3*ncv*ncv+5*ncv;   //size of workl array 
   int *iparam = (int *) malloc(11*sizeof(int));
   int *ipntr  = (int *) malloc(14*sizeof(int));
   int *select = (int *) malloc(ncv*sizeof(int));
   spinor *workd;
   _Complex double *workev, *resid, *workl; 
   double *rwork,*rd;
   workd  = (spinor *) alloc_aligned_mem(3*ldv*sizeof(spinor));
   workev = (_Complex double *) alloc_aligned_mem(3*ncv*sizeof(_Complex double));
   resid  = (_Complex double *) alloc_aligned_mem(12*n*sizeof(_Complex double));
   workl  = (_Complex double *) alloc_aligned_mem(lworkl*sizeof(_Complex double));
   rwork  = (double *) alloc_aligned_mem(  ncv*sizeof(double));
   rd     = (double *) alloc_aligned_mem(3*ncv*sizeof(double));

   char bmat[]= "I";     /* Specifies that the right hand side matrix
                            should be the identity matrix; this makes
                            the problem a standard eigenvalue problem.
                         */
   char cA='A';

   int parallel;
   #ifdef MPI
   MPI_Comm comm; //communicator used when we call PARPACK
   int comm_err = MPI_Comm_dup(MPI_COMM_WORLD,&comm); //duplicate the MPI_COMM_WORLD to create a communicator to be used with arpack
   if(comm_err != MPI_SUCCESS) { //error when trying to duplicate the communicator
     if(g_proc_id == g_stdio_proc){
       fprintf(stderr,"MPI_Comm_dup return with an error. Exciting...\n");
       return -1;
     }
   }
   #endif


   #ifdef MPI
     parallel=1;
   #else
     parallel=0;
   #endif
   
   //set params for arpack
   ido     = 0;                 //reverse communication parameter that intially is set to zero
   (*info) = 0;                 //means use a random starting vector with Arnoldi
   ishfts  = 1;                 //use exact shifts (other options can be found in the documentation)
   mode    = 1;                 //regular mode
   iparam[0] = ishfts;          
   iparam[2] = maxiter;         
   iparam[6] = mode;             

   /*
     M A I N   L O O P (Reverse communication)  
   */
   do
   {
      #ifndef MPI (serial code)
      _FT(znaupd)(&ido, bmat, &N, which, &nev, &tol, resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, (_Complex double *) workd, 
                  workl, &lworkl,rwork,info );
      #else
      _FT(pznaupd)(&comm, &ido, bmat, &N, which, &nev, &tol, resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, (_Complex double *) workd, 
                  workl, &lworkl,rwork,info );
      #endif
      if ((ido==-1)||(ido==1)){ 
         av(workd+(ipntr[1]-1)/12, workd+(ipntr[0]-1)/12);
      }
      
   } while ((ido==-1)||(ido==1));
   
/*
 Check for convergence 
*/
     if ( (*info) < 0 ) 
     {
         if(g_proc_id == g_stdio_proc){
            fprintf(stderr,"Error with _naupd, info = %d\n", *info);
            fprintf(stderr,"Check the documentation of _naupd\n");}
     }
     else 
     { 
        //compute eigenvectors if desired
        #ifndef MPI
        _FT(zneupd) (&comp_evecs,&cA, select,evals,v,&N,&sigma, 
                     workev,bmat,&N,which,&nev,&tol,resid,&ncv, 
                     v,&N,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,&ierr);
        #else
        _FT(pzneupd) (&comm,&comp_evecs,&cA,select,evals,v,&N,&sigma, 
                      workev,bmat,&N,which,&nev,&tol,resid,&ncv, 
                      v,&N,iparam,ipntr,workd,workl,&lworkl, 
                      rwork,&ierr);
        #endif

        /*
        %----------------------------------------------%
        | Eigenvalues are returned in the one          |
        | dimensional array evals.  The corresponding  |
        | eigenvectors are returned in the first NCONV |
        | (=IPARAM[4]) columns of the two dimensional  | 
        | array V if requested.  Otherwise, an         |
        | orthogonal basis for the invariant subspace  |
        | corresponding to the eigenvalues in evals is |
        | returned in V.                               |
        %----------------------------------------------%
        */

        if( ierr!=0) 
        {
           if(g_proc_id == g_stdio_proc){
             fprintf(stderr,"Error with _neupd, info = %d \n",ierr);
             fprintf(stderr,"Check the documentation of _neupd. \n");}
        }
        else //report eiegnvalues and their residuals
        {
             (*nconv) = iparam[4];
             for(j=0; j< (*nconv); j++)
             {
               /*
               %---------------------------%
               | Compute the residual norm |
               |                           |
               | ||A*x - lambda*x ||/||x|| |
               |                           |
               | for the NCONV accurately  |
               | computed eigenvalues and  |
               | eigenvectors.  (iparam[4] |
               | indicates how many are    |
               | accurate to the requested |
               | tolerance)                |
               %---------------------------%
               */
               
               //compute the residual
               //IMPORTANT: our eigenvectors are of lattice size. In order to apply the Dirac operator
               //we need to copy them to a spinor with lattice size + buffer size needed for communications
               assign(workd+ldv,&v[j*n],n);
               av(workd,workd+ldv);
               assign_diff_mul(workd,&v[j*n], evals[j], n);
               rd[j*3]   = creal(evals[j]);
               rd[j*3+1] = cimag(evals[j]);
               rd[j*3+2] = sqrt(square_norm(workd,n,parallel));
               rd[j*3+2] = rd[j*3+2] /sqrt(square_norm(&v[j*n],n,parallel));

               if(g_proc_id == g_stdio_proc)
                  fprintf(stdout,"RitzValue %d  %f  %f  ||A*x-lambda*x||/||x||= %f\n",j,rd[j*3],rd[j*3+1],rd[j*3+2]);
              }
        }

        /*Print additional convergence information.*/
        if( (*info)==1)
        {
           if(g_proc_id == g_stdio_proc)
             fprintf(stderr,"Maximum number of iterations reached.\n");
        }
        else
        { 
           
           if(g_proc_id == g_stdio_proc)
           {
              if((*info)==3)
              {  
                 fprintf(stderr,"No shifts could be applied during implicit\n");
                 fprintf(stderr,"Arnoldi update, try increasing NCV.\n");
              }
         
              fprintf(stdout,"_NDRV1\n");
              fprintf(stdout,"=======\n");
              fprintf(stdout,"Size of the matrix is %d\n", 12*n);
              fprintf(stdout,"The number of Ritz values requested is %d\n", nev);
              fprintf(stdout,"The number of Arnoldi vectors generated is %d\n", ncv);
              fprintf(stdout,"What portion of the spectrum: %s\n", which);
              fprintf(stdout,"The number of converged Ritz values is %d\n", (*nconv) ); 
              fprintf(stdout,"The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
              fprintf(stdout,"The number of OP*x is %d\n", iparam[8]);
              fprintf(stdout,"The convergence criterion is %f\n", tol);
           }
        }
     }  //if(info < 0) else part
     
     return;
}

void evals_arpack_poly_hermitian(int n, int nev, int ncv, char *which, _Complex double *evals, spinor *v, double tol, int maxiter, 
                                 matrix_mult av, int cheb_k, double evmin, double evmax, int *info, int *nconv)
{
   int N,ldv;
  
   if(n==VOLUME) //full 
     ldv= VOLUMEPLUSRAND;
   else          //even-odd
     ldv= VOLUMEPLUSRAND/2;
   
   N=12*n;  //size of the matrix as complex double

   //check input
   if(nev>=N) nev=N-1;
   if(ncv < (nev+1)) ncv = nev+1;

 
   //parameters for the ARPACK routines--see documentation or source files
   int ido,lworkl,j,ierr,ishfts,mode;
   _Complex double sigma;
   int comp_evecs=1; //option to compute the eigenvectors (always set to true here)
   lworkl  = 3*ncv*ncv+5*ncv;   //size of workl array 
   int *iparam = (int *) malloc(11*sizeof(int));
   int *ipntr  = (int *) malloc(14*sizeof(int));
   int *select = (int *) malloc(ncv*sizeof(int));
   spinor *workd;
   _Complex double *workev, *resid, *workl; 
   double *rwork,*rd;
   workd  = (spinor *) alloc_aligned_mem(3*ldv*sizeof(spinor));
   workev = (_Complex double *) alloc_aligned_mem(3*ncv*sizeof(_Complex double));
   resid  = (_Complex double *) alloc_aligned_mem(12*n*sizeof(_Complex double));
   workl  = (_Complex double *) alloc_aligned_mem(lworkl*sizeof(_Complex double));
   rwork  = (double *) alloc_aligned_mem(  ncv*sizeof(double));
   rd     = (double *) alloc_aligned_mem(3*ncv*sizeof(double));

   char bmat[]= "I";     /* Specifies that the right hand side matrix
                            should be the identity matrix; this makes
                            the problem a standard eigenvalue problem.
                         */
   char cA='A';

   int parallel;
   #ifdef MPI
   MPI_Comm comm; //communicator used when we call PARPACK
   int comm_err = MPI_Comm_dup(MPI_COMM_WORLD,&comm); //duplicate the MPI_COMM_WORLD to create a communicator to be used with arpack
   if(comm_err != MPI_SUCCESS) { //error when trying to duplicate the communicator
     if(g_proc_id == g_stdio_proc){
       fprintf(stderr,"MPI_Comm_dup return with an error. Exciting...\n");
       return -1;
     }
   }
   #endif


   #ifdef MPI
     parallel=1;
   #else
     parallel=0;
   #endif
   
   //set params for arpack
   ido     = 0;                 //reverse communication parameter that intially is set to zero
   (*info) = 0;                 //means use a random starting vector with Arnoldi
   ishfts  = 1;                 //use exact shifts (other options can be found in the documentation)
   mode    = 1;                 //regular mode
   iparam[0] = ishfts;          
   iparam[2] = maxiter;         
   iparam[6] = mode;             

   /*
     M A I N   L O O P (Reverse communication)  
   */
   do
   {
      #ifndef MPI (serial code)
      _FT(znaupd)(&ido, bmat, &N, which, &nev, &tol, resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, (_Complex double *) workd, 
                  workl, &lworkl,rwork,info );
      #else
      _FT(pznaupd)(&comm, &ido, bmat, &N, which, &nev, &tol, resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, (_Complex double *) workd, 
                  workl, &lworkl,rwork,info );
      #endif
      if ((ido==-1)||(ido==1)){ 
         //av(workd+(ipntr[1]-1)/12, workd+(ipntr[0]-1)/12);
         cheb_poly_precon_op(workd+(ipntr[1]-1)/12, workd+(ipntr[0]-1)/12,av,n,evmin,evmax,cheb_k);
      }
      
   } while ((ido==-1)||(ido==1));
   
/*
 Check for convergence 
*/
     if ( (*info) < 0 ) 
     {
         if(g_proc_id == g_stdio_proc){
            fprintf(stderr,"Error with _naupd, info = %d\n", *info);
            fprintf(stderr,"Check the documentation of _naupd\n");}
     }
     else 
     { 
        //compute eigenvectors if desired
        #ifndef MPI
        _FT(zneupd) (&comp_evecs,&cA, select,evals,v,&N,&sigma, 
                     workev,bmat,&N,which,&nev,&tol,resid,&ncv, 
                     v,&N,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,&ierr);
        #else
        _FT(pzneupd) (&comm,&comp_evecs,&cA,select,evals,v,&N,&sigma, 
                      workev,bmat,&N,which,&nev,&tol,resid,&ncv, 
                      v,&N,iparam,ipntr,workd,workl,&lworkl, 
                      rwork,&ierr);
        #endif

        /*
        %----------------------------------------------%
        | Eigenvalues are returned in the one          |
        | dimensional array evals.  The corresponding  |
        | eigenvectors are returned in the first NCONV |
        | (=IPARAM[4]) columns of the two dimensional  | 
        | array V if requested.  Otherwise, an         |
        | orthogonal basis for the invariant subspace  |
        | corresponding to the eigenvalues in evals is |
        | returned in V.                               |
        %----------------------------------------------%
        */

        if( ierr!=0) 
        {
           if(g_proc_id == g_stdio_proc){
             fprintf(stderr,"Error with _neupd, info = %d \n",ierr);
             fprintf(stderr,"Check the documentation of _neupd. \n");}
        }
        else //report eiegnvalues and their residuals
        {
             (*nconv) = iparam[4];
             for(j=0; j< (*nconv); j++)
             {
               /*
               %---------------------------%
               | Compute the residual norm |
               |                           |
               | ||A*x - lambda*x ||/||x|| |
               |                           |
               | for the NCONV accurately  |
               | computed eigenvalues and  |
               | eigenvectors.  (iparam[4] |
               | indicates how many are    |
               | accurate to the requested |
               | tolerance)                |
               %---------------------------%
               */
               
               //compute the residual
               //IMPORTANT: our eigenvectors are of lattice size. In order to apply the Dirac operator
               //we need to copy them to a spinor with lattice size + buffer size needed for communications
               assign(workd+ldv,&v[j*n],n);
               //av(workd,workd+ldv);
               cheb_poly_precon_op(workd,workd+ldv,av,n,evmin,evmax,cheb_k);
               assign_diff_mul(workd,&v[j*n], evals[j], n);
               rd[j*3]   = creal(evals[j]);
               rd[j*3+1] = cimag(evals[j]);
               rd[j*3+2] = sqrt(square_norm(workd,n,parallel));
               rd[j*3+2] = rd[j*3+2] /sqrt(square_norm(&v[j*n],n,parallel));
               if(g_proc_id == g_stdio_proc)
                  fprintf(stdout,"RitzValue %d  %f  %f  ||A*x-lambda*x||/||x||= %f\n",j,rd[j*3],rd[j*3+1],rd[j*3+2]);
              }
        }

        /*Print additional convergence information.*/
        if( (*info)==1)
        {
           if(g_proc_id == g_stdio_proc)
             fprintf(stderr,"Maximum number of iterations reached.\n");
        }
        else
        { 
           
           if(g_proc_id == g_stdio_proc)
           {
              if((*info)==3)
              {  
                 fprintf(stderr,"No shifts could be applied during implicit\n");
                 fprintf(stderr,"Arnoldi update, try increasing NCV.\n");
              }
         
              fprintf(stdout,"_NDRV1\n");
              fprintf(stdout,"=======\n");
              fprintf(stdout,"Size of the matrix is %d\n", 12*n);
              fprintf(stdout,"The number of Ritz values requested is %d\n", nev);
              fprintf(stdout,"The number of Arnoldi vectors generated is %d\n", ncv);
              fprintf(stdout,"What portion of the spectrum: %s\n", which);
              fprintf(stdout,"The number of converged Ritz values is %d\n", (*nconv) ); 
              fprintf(stdout,"The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
              fprintf(stdout,"The number of OP*x is %d\n", iparam[8]);
              fprintf(stdout,"The convergence criterion is %f\n", tol);
              fprintf(stdout,"The system is preconditioned using chebyshev polynomils of degree %d\n",cheb_k);
              fprintf(stdout,"Estimates of the eigenvalue bounds were given as [ %f , %f]\n",evmin,evmax); 
           }
        }
     }  //if(info < 0) else part
     
     return;
}



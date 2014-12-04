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
  double tol, 
  int maxiter, 
  matrix_mult av, 
  int *info, 
  int *nconv)
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
  cheb_k   : (IN) degree of the chebyshev polynomial to be used for acceleration (irrelevant when use_acc=0 and init_resid_arpack=0)
  amin,amax: (IN) bounds of the interval [amin,amax] for the acceleration polynomial (irrelevant when use_acc=0 and init_resid_arpack=0)
  evals : (OUT) array of size nev+1 which has the computed nev Ritz values
  v     : computed eigenvectors. Size is n*ncv spinors.
  tol    : Requested tolerance for the accuracy of the computed eigenvectors.
           A value of 0 means machine precision.
  maxiter: maximum number of restarts (iterations) allowed to be used by ARPACK
  av     : operator for computing the action of the matrix on the vector
           av(vout,vin) where vout is output spinor and vin is input spinors.
  info   : output from arpack. 0 means that it converged to the desired tolerance. 
           otherwise, an error message is printed to stderr 
  nconv  : actual number of converged eigenvectors.
*/ 
{
   //create the MPI communicator
   #ifdef MPI
   MPI_Comm comm; //communicator used when we call PARPACK
   int comm_err = MPI_Comm_dup(MPI_COMM_WORLD,&comm); //duplicate the MPI_COMM_WORLD to create a communicator to be used with arpack
   if(comm_err != MPI_SUCCESS) { //error when trying to duplicate the communicator
     if(g_proc_id == g_stdio_proc){
       fprintf(stderr,"MPI_Comm_dup return with an error. Exciting...\n");
       exit(-1);
     }
   }
   #endif

   int parallel;
   #ifdef MPI
     parallel=1;
   #else
     parallel=0;
   #endif

   int ido=0;           //control of the action taken by reverse communications
                        //set initially to zero

   char bmat= 'I';     /* Specifies that the right hand side matrix
                          should be the identity matrix; this makes
                          the problem a standard eigenvalue problem.
                       */

   //matrix dimensions 
   int ldv,N,LDV;
  
   if(n==VOLUME) //full 
     ldv= VOLUMEPLUSRAND;
   else          //even-odd
     ldv= VOLUMEPLUSRAND/2;

   //dimesnions as complex variables
   N  =12*n;       //dimension
   LDV=12*ldv;   //leading dimension (including communication buffers)
   


   char *which_evals;

   if(which==0)
     which_evals="SR";
   if(which==1)
     which_evals="LR";
   if(which==2)
     which_evals="SM";
   if(which==3)
     which_evals="LM";


   //check input
   if(nev>=N) nev=N-1;
   if(ncv < (nev+1)) ncv = nev+1;

   spinor *resid  = (spinor *) alloc_aligned_mem(ldv*sizeof(spinor));

   int *iparam = (int *) malloc(11*sizeof(int));

   if((use_acc==0) || (use_acc==1))
     iparam[0]=1;
   if(use_acc==2)
     iparam[0]=0;

   iparam[2]=maxiter;

   iparam[3]=1;

   iparam[6]=1;


   int *ipntr  = (int *) malloc(14*sizeof(int));

   _Complex double *workd  = (_Complex double *) alloc_aligned_mem(3*N*sizeof(_Complex double)); 

   int lworkl=3*ncv*ncv+5*ncv;

   _Complex double *workl=(_Complex double *) alloc_aligned_mem(lworkl*sizeof(_Complex double));

   double *rwork  = (double *) alloc_aligned_mem(  ncv*sizeof(double));

   int rvec=1; //always compute eigenvectors

   char howmany='A';   //compute eigenvectors

   //spinor *zv; //this is for the eigenvectors and won't be referenced when howmany='P'

   int *select = (int *) malloc(ncv*sizeof(int)); //since all Ritz vectors or Schur vectors are computed no need to initialize this array

   _Complex double sigma;
    
   _Complex double *workev = (_Complex double *) alloc_aligned_mem(2*ncv*sizeof(_Complex double));
   
   double d1,d2,d3;

   //if(init_resid_arpack==0)
      (*info) = 0;                 //means use a random starting vector with Arnoldi

   spinor *vin   = (spinor *) alloc_aligned_mem(ldv*sizeof(spinor)); //input spinor 
   spinor *vout  = (spinor *) alloc_aligned_mem(ldv*sizeof(spinor)); //output spinor

   int i,j;



   /*
     M A I N   L O O P (Reverse communication)  
   */
   do
   {
      #ifndef MPI 
      _FT(znaupd)(&ido, &bmat, &N, which_evals, &nev, &tol, (_Complex double *) resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, workd, 
                  workl, &lworkl,rwork,info );
      #else
      _FT(pznaupd)(&comm, &ido, &bmat, &N, which_evals, &nev, &tol, (_Complex double *) resid, &ncv,
                  (_Complex double *) v, &N, iparam, ipntr, workd, 
                  workl, &lworkl,rwork,info );
      #endif
      if ((ido==-1)||(ido==1)){

         assign(vin, (spinor *) workd+(ipntr[0]-1)/12,n);

         if((use_acc==0) || (use_acc==2) )
           av(vout,vin);
         else 
           cheb_poly_op(vout,vin,av,n,amin,amax,cheb_k);

         assign((spinor *) workd+(ipntr[1]-1)/12, vout, n);
      }

      if( (ido==3) & (iparam[0]==0)){ //provide implicit shifts as roots of the chebyshev polynoimial
         cheb_poly_roots(&workl[ipntr[13]],iparam[7],amin,amax);
      }


      
   } while (ido != 99);
   
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
        (*nconv) = iparam[4];
        if(g_proc_id == g_stdio_proc){
          fprintf(stderr,"number of converged eigenvectors = %d\n", *nconv);}

        //compute eigenvectors 
        #ifndef MPI
        _FT(zneupd) (&rvec,&howmany, select,evals,(_Complex double *) v,&N,&sigma, 
                     workev,&bmat,&N,which_evals,&nev,&tol,(_Complex double *) resid,&ncv, 
                     (_Complex double *) v,&N,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,info);
        #else
        _FT(pzneupd) (&comm,&rvec,&howmany, select,evals, (_Complex double *) v,&N,&sigma, 
                     workev,&bmat,&N,which_evals,&nev,&tol,(_Complex double *) resid,&ncv, 
                     (_Complex double *) v,&N,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,info);
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

        if( (*info)!=0) 
        {
           if(g_proc_id == g_stdio_proc){
             fprintf(stderr,"Error with _neupd, info = %d \n",(*info));
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
               assign(vin,&v[j*n],n);
               if((use_acc==0) || (use_acc==2))
                 av(vout,vin);
               else
                 cheb_poly_op(vout,vin,av,n,amin,amax,cheb_k);

               assign_diff_mul(vout,&v[j*n], evals[j], n);
               
               d1 = sqrt(square_norm(vout,n,parallel));
               d2 = sqrt(square_norm(&v[j*n],n,parallel));

               if(g_proc_id == g_stdio_proc)
                  fprintf(stdout,"RitzValue %d  %+e  %+e  error= %+e \n",j,creal(evals[j]),cimag(evals[j]),d1/d2);
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
              fprintf(stdout,"Size of the matrix is %d\n", N);
              fprintf(stdout,"The number of Ritz values requested is %d\n", nev);
              fprintf(stdout,"The number of Arnoldi vectors generated is %d\n", ncv);
              fprintf(stdout,"What portion of the spectrum: %s\n", which_evals);
              fprintf(stdout,"The number of converged Ritz values is %d\n", (*nconv) ); 
              fprintf(stdout,"The number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
              fprintf(stdout,"The number of OP*x is %d\n", iparam[8]);
              fprintf(stdout,"The convergence criterion is %f\n", tol);
           }
          
        }
     }  //if(info < 0) else part
     
     return;
}


/*********************************************************************** 
 * preconditioning related functions.
 * Note: there exist already polynomial preconditioning operators
 *       related to the chebychev polynomials in poly_precon.h files
 *       that was developed earlier by Carsten
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/

#include "solver/precon.h"

void cheb_poly_precon(spinor * const R, spinor * const S, matrix_mult f, const int N, const double evmin, const double evmax, const int k)
{

/*
 R = T_{k+1}(Q)*S (k>=-1) where S is input spinor and Q is the operator given by the matrix-vector multipliction where Q*v is given by the
 matrix-vector multiplication opertor f. T_{k+1}(Q) is the chebychev polynomial given by:
  
       evmin and evmax are minimum and maximum eigenvalues of the operator under consideration (assumed to be Hermitian here)  

       theta = (evmax+evmin)/2,   delta=(evmax-evmin)/2
       sigma_0 =1 ,     sigma_1=theta/delta
       rho_{-1}=0 ,     rho_0=1/(2*sigma_1)
       rho_k   = 1/(2*sigma_1 - rho_{k-1})

       T_{-1}(t) = 0,    T_0(t)=1,  T_1(t)=1-t/theta
       
       T_{k+1}(t) = rho_k*{ 2*(sigma_1-t/delta)*T_k(t) - rho_{k-1}*T_{k-1}) 

R: output spinor
S: input spinor
f: matrix-vector multiplication operator
N: size of the spinor
evmin: minimum value of the interval (approximate lowest eigenvalue)
evmax: maximum value of the interval (approximate largest eigenvalue)
k: order of the Chebeychev polynomial is k+1
*/

   double theta,delta,sigma1,rho_prev,rho;
   double d1,d2,d3;
   static int initp=0;
   static spinor *vm1,*tmpv1,*tmpv2;

   if(k < -1){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest order of the polynomial is 0 (k=-1).\n");
        exit(1);
    }

    //T_0(Q)=1 
    assign(R,S,N);
    if(k== -1){
      return;
    }


    theta=0.5*(evmax+evmin);
    delta=0.5*(evmax-evmin);

    if(delta < 5e-15){ //check the size of the interval
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: evmax-evmin < 5e-15.\n");
        exit(1);
    }

    sigma1=theta/delta;
    rho_prev=0.0;

    int LDN;
    if(N==VOLUME)
       LDN = VOLUMEPLUSRAND;
    else
       LDN = VOLUMEPLUSRAND/2;

    if(initp==0)
    {
       //alloacate memory for needed spinors
       vm1   = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       tmpv1 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       tmpv2 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       
       initp=1;
    }




    zero_spinor_field(vm1,N);

    for(int i=0; i <= k; i++)
    {
       //rho_{-1}=0 , rho_0=1/(2*sigma_1)
       //rho_k   = 1/(2*sigma_1 - rho_{k-1})
       //T_{i+1}(t) = rho_i*{ 2*(sigma_1-t/delta)*T_i(t) - rho_{i-1}*T_{i-1}) 

       rho=1.0/(2*sigma1-rho_prev);

       d1=  2*rho*sigma1;
       d2= -2.0*rho/delta;
       d3= -rho*rho_prev;

       f(tmpv1,R);
       assign(tmpv2,R,N);

       assign_mul_add_mul_add_mul_r(R,tmpv1,vm1,d1,d2,d3,N);
       assign(vm1,tmpv2,N);

       rho_prev=rho;
    }

    return;

}


void cheb_poly_precon_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double evmin, const double evmax, const int k)
{

   static int initop=0;
   static spinor *v;
   int LDN;
   if(N==VOLUME)
       LDN = VOLUMEPLUSRAND;
    else
       LDN = VOLUMEPLUSRAND/2;

   if(initop==0)
   {
      v = (spinor *) alloc_ligned_mem(LDN*sizeof(spinor));
      initop=1;
   }

   f(v,S);
   cheb_poly_precon(R,v,f,N,evmin,evmax,k);
   return;
}

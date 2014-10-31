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
 R = T_{k}(Q)*S (k>=0) where S is input spinor and Q is the operator given by the matrix-vector multipliction where Q*v is given by the
 matrix-vector multiplication opertor f. T_{k}(Q) is the chebychev polynomial.
  
 For -1 =< t =< 1, the Chebyshev polynomials are given by

       T_0(t) = 1    
       T_1(t) = t
       T_j(t) = 2*t*T_{j-1}(t) - T_{j-2}(t) where j=2,3,...

and the roots of T_n(t) in the interval [-1,1] for n>=1 are given by

       t_i = cos(pi/2 *(2*i-1)/j) for i=1,2,..,j   and j>=1

In our case, we choose an interval [evmin,evmax] which will require shifting the polynomial.

So we define t = 2/(evmax-evmin)*x - (evmax+evmin)/(evmax-evmin)

where x here is our matrix Q.

So, as a polynomial of Q, we have 
    T_0(Q) = 1

    T_1(Q) = 2/(evmax-evmin)*Q - (evmax+evmin)/(evmax-evmin)

    T_j(Q) = 4/(evmax-evmin)*Q*T_{j-1}(Q) - 2*(evmax+evmin)/(evmax-evmin)*T_{j-1}(Q) - T_{j-2}(Q)  


and the roots of this polynomial are 

    q_i = (evmax+evmin)/2+ (evmax-evmin)/2*cos(pi/2 *(2*i-1)/j)   for i=1,2,..,j and j>=1 


R: output spinor
S: input spinor
f: matrix-vector multiplication operator
N: size of the spinor
evmin: minimum value of the interval
evmax: maximum value of the interval
k: order of the Chebeychev polynomial
*/

   double d1,d2,d3;
   static int initp=0;
   static spinor *tmpv1,*tmpv2,*tmpv3;

   if(k < 0){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest allowed order of the polynomial is 0.\n");
        exit(1);
   }


   //T_0(Q)=1 
   assign(R,S,N);
   if(k== 0){
      return;
   }


   //T_1(Q) = 2/(evmax-evmin)*Q - (evmax+evmin)/(evmax-evmin)
   d1 = 2.0/(evmax-evmin);
   d2 = -(evmax+evmin)/(evmax-evmin);
   f(R,S); //R=Q(S)
   assign_mul_add_mul_r(R,S,d1,d2,N);
   if(k==1){
     return;
   }
   
   //degree >=2
   //==========
   int LDN;
   if(N==VOLUME)
      LDN = VOLUMEPLUSRAND;
   else
      LDN = VOLUMEPLUSRAND/2;

   //allocate needed memory
   if(initp==0)
   {
       tmpv1 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       tmpv2 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       tmpv3 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       initp=1;
   }


   //T_j(Q) = 4/(evmax-evmin)*Q*T_{j-1}(Q) - 2*(evmax+evmin)/(evmax-evmin)*T_{j-1}(Q) - T_{j-2}(Q)  
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   d1 = 4.0/(evmax-evmin);
   d2 = -2.0*(evmax+evmin)/(evmax-evmin);
   d3 = -1.0;

   

    assign(tmpv1,S,N);
    assign(tmpv2,R,N);
    for(int i=2; i <= k; i++)
    {
       f(R,tmpv2);
       assign_mul_add_mul_add_mul_r(R,tmpv2,tmpv1,d1,d2,d3,N);
       assign(tmpv1,tmpv2,N);
       assign(tmpv2,R,N);
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
      v = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
      initop=1;
   }

   f(v,S);
   cheb_poly_precon(R,v,f,N,evmin,evmax,k);
   return;
}

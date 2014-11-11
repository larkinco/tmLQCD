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


//compute R= T_k(A) S where T_k(A) = C_k(d(A))/C_k(d(0)) and d(A) = 2A/(b-a)-(b+a)/(b-a) where [a,b] are the boundaries of a given interval
//R: output spinor (out)
//S: input spinor (in) 
//f: matrix-vector multiplication function corresponding to the operator A (in)
//N: size of the spinors      (in)
//a,b: limits of the interval (in)
//k: degree of the polynomial requested >=1 (in)
void cheb_poly_precon_residual(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k)
{
/*
This is the polynomial T_{k}(x) = 1 -xP_{k-1}(x) defined such that |1-xP_{k-1}(x)| is minimized over the interval [a,b]. P_{k-1}(x) is the approximate
inverse over this interval and T(x) is the residual polynomial. When solving a linear system Ay=b with initial residual r_0=b-Ay_0, we see that the 
resiudal of the preconditioned system P(A)A y = P(A)b has a residual T(A)r_0. For eigenvalue problems, we can use T_{k}(A) as a filter to minimize
the contribution coming from eigenvalues in the interval [a,b] and enhance the other part of the spectrum. This will be used with ARPACK for example.
For the linear system on the other hand we will use P(A) as a preconditioner and it is defined in the function cheb_poly_precon_op. Note that the roots
of T(x) are related to the roots of the Chebyshev polynomial of the first kind.

Propertires of the Chebyshev polynomials of the first kind C_k(x) for x in the interval [-1,1] and k=0,1,2,3,...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

k=0: C_0(x) =1
k=1: C_1(x) = x
...
three-term recurrence:   C_{k+1}(x) = 2xC_k(x) - C_{k-1}

C_k(x) has k simple roots in the interval [-1,1] given by
x_l = cos(pi/2 (2l-1)/k) for k=1,2,3,... and l=1,2,..,k

In the interval [a,b] we can define a shifted polynomial by C_k(d(x) where d(x) = 2x/(b-a) - (b+a)/(b-a)
such that d(a)=-1, and d(b)=+1. 

d(x) = x/delta - theta/delta  with delta = (b-a)/2 and theta=(b+a)/2

R_k(x) = C_k(d(x) is the shifted Chebyshev polynomial with 

R_0(x) =1
R_1(x) = d(x)

R_{k+1}(x) = 2d(x)R_k(x)-R_{k-1}(x)

T_k(x) is the normalized shifted Chebyshev polynomial given by:

For k=1,2,3,....   T_k(x) = R_k(x)/sigma_k(x)

with sigma_0 = 1, sigma_1=d(0) = - theta/delta
sigma_{k+1} =  -2theta/delta sigma_k - sigma_{k-1}

These relations define the polynomial T_k(x)

The roots of T_k(x) are the same as the roots of C_k(d(x)) and are given by

x_l = (b-a)/2*[cos(pi/2 (2*l-1)/k)+(b+a)/(b-a)]

*/
   static int initp=0;
   static spinor *v1,*v2,*v3;

   double delta,theta;
   double sigma,sigma1,sigma2;
   double d1,d2,d3;

   if(k < 0){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest allowed order of the residual polynomial is 0.\n");
        exit(1);
   }

   delta = (b-a)/2.0;
   theta = (b+a)/2.0;

   sigma  = 1.0;
   sigma1 = -theta/delta;

   //T_0(Q)=1 
   assign(R,S,N);
   if(k== 0){
      return;
   }


   //T_1(Q) = [2/(b-a)*Q - (b+a)/(b-a)] / sigma_1 = -Q/theta +1 
   d1 = -1.0/theta;
   d2 =  1.0;
   f(R,S); //R=Q(S)
   assign_mul_add_mul_r(R,S,d1,d2,N);
   if(k==1){
     return;
   }
   
   //degree >=2
   //==========

   //T_0 = S
   //T_1 = R


   int LDN;
   if(N==VOLUME)
      LDN = VOLUMEPLUSRAND;
   else
      LDN = VOLUMEPLUSRAND/2;

   //allocate needed memory
   if(initp==0)
   {
       v1 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       v2 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       v3 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       initp=1;
   }


   // T_j(Q) = 1/(sigma_j*delta) Q - theta/(delta*sigma_j) T_{j-1}(Q) -1/sigma_j T_{j-2}(Q)  for j=2,3,...  
   //-----------------------------------------------------------------------------------------------------

   assign(v1,S,N);
   assign(v2,R,N);

   //v1 = T_0,  v2 = T_1


   for(int i=2; i <= k; i++)
   {
      sigma2 = -2.0*theta/delta*sigma1 - sigma;
      
      d1 = 1.0/(sigma2*delta);
      d2 = -theta/(sigma2*delta);
      d3 = -1.0/sigma2;

      f(R,v2);
      assign_mul_add_mul_add_mul_r(R,v2,v1,d1,d2,d3,N);

      assign(v1,v2,N);
      assign(v2,R,N);

      sigma  = sigma1;
      sigma1 = sigma2;  
   }

   return;

}


//This is the version of the preconditioner to be used with linear system solution
void cheb_poly_precon_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k)
/*
Computes the action of the polynomial P_k(x) which minimizes the |1-xP_k(x)| over the interval [a,b]. This polynomial is used as 
an approximation to the inverse of the matrix over that interval and can be used as a preconditioner for CG. The recurrence relations
could be derived similar to the residual polynomial above and are given by:

sigma_0 = 1, sigma_1=d(0) = - theta/delta
sigma_{k+1} =  -2theta/delta sigma_k - sigma_{k-1}

P_0(Q) = -1/theta

P_1(Q) = - (4*theta+2*Q)/(2*theta^2-delta^2)

P_{k+1}(Q) = -(2*sigma_{k+1}/(delta*sigma_{k+2}))- 2*sigma_{k+1}/(delta*sigma_{k+2})*(theta-Q)*P_k(Q) - sigma_k/sigma_{k+2}*P_{k-1}(Q)

*/
{
   static int initp=0;
   static spinor *v1,*v2,*v3,*v4;

   double delta,theta;
   double sigma0,sigma1,sigma2,sigma3;
   double d1,d2,d3,d4;

   if(k < 0){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest allowed order of the residual polynomial is 0.\n");
        exit(1);
   }

   int LDN;
   if(N==VOLUME)
      LDN = VOLUMEPLUSRAND;
   else
      LDN = VOLUMEPLUSRAND/2;

   //allocate needed memory
   if(initp==0)
   {
       v1 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       v2 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       v3 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       v4 = (spinor *) alloc_aligned_mem(LDN*sizeof(spinor));
       initp=1;
   }

   delta = (b-a)/2.0;
   theta = (b+a)/2.0;

   sigma0 =  1.0;
   sigma1 = -theta/delta;
   sigma2 = -2.0*theta/delta*sigma1 - sigma0;

   //P_0(Q) = -1/theta;
   d1 = -1.0/theta;
   mul_r(v1,d1,S,N);
   if(k== 0){
      return;
   }


   //P_1(Q) = -(4*theta+2*Q)/(2*theta^2-delta^2) 
   d1 = -4.0*theta/(2.0*theta*theta-delta*delta);
   d2 = -2.0/(2*theta*theta-delta*delta);
   f(R,S); //R=Q(S)
   assign_mul_add_mul_r(R,S,d1,d2,N);
   if(k==1){
     return;
   }
   
   //P_{k+1}(Q) = -(2*sigma_{k+1}/(delta*sigma_{k+2}))- 2*sigma_{k+1}/(delta*sigma_{k+2})*(theta-Q)*P_k(Q) - sigma_k/sigma_{k+2}*P_{k-1}(Q)
   //--------------------------------------------------------------------------------------------------------------------------------------
   assign(v2,R,N);

   //v1 = P_0,  v2 = P_1
   for(int i=2; i <= k; i++)
   {
      sigma3 = -2.0*theta/delta*sigma2 - sigma1;
      
      d1 = +2.0*sigma2/delta/sigma3;
      d2 = -2.0*sigma2*theta/delta/sigma3;
      d3 = -d2;
      d4 = -sigma1/sigma3;

      f(R,v2);
      assign_mul_add_mul_add_mul_add_mul_r(R,v2,S,v1,d1,d2,d3,d4,N);

      
      sigma0 = sigma1;
      sigma1 = sigma2;
      sigma2 = sigma3;

      assign(v1,v2,N);
      assign(v2,R,N);  
   }

   return;
}

void cheb_poly_roots(_Complex double *roots, const int k, const double a, const double b)
/*
roots of the shifted Chebyshev polynomial in the interval [a,b]
The roots of C_k(d(x)) and are given by
x_l = (b-a)/2*[cos(pi/2 (2*l-1)/k)+(b+a)/(b-a)]
*/
{
   double PI=3.141592653589793;

   double d1,d2,d3;

   d1=0.5*(b+a);
   d2=0.5*(b-a);
   d3=PI/(double)k;


   int i;

   for(i=1; i<=k; i++)
      roots[i-1] = d1+d2*cos(d3*(i-0.5));

   
   return ;
}






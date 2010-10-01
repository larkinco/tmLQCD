/***********************************************************************
 *
 * Copyright (C) 2010 Florian Burger
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
 *  
 * File: mixed_solve.cu
 *
 * CUDA GPU mixed_solver for EO and non-EO
 * CUDA kernels for Hopping-Matrix and D_tm
 *
 * The externally accessible functions are
 *
 *
 *   extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec, const int N)
 *
 *  extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
           double eps, const int rel_prec,const int N)
 *
 * input:
 *   Q: source
 * inout:
 *   P: initial guess and result
 * 
 *
 **************************************************************************/



#include <cuda.h>
#include <cuda_runtime.h>
#include "cublas.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "../global.h"
#include "cudaglobal.h"
#include "mixed_solve.h"
#include "cudadefs.h"
#include <math.h>


extern "C" {
#include "../tm_operators.h"
#include "../linalg_eo.h"
#include "../start.h"
#include "../complex.h"
#include "../read_input.h"
#include "../geometry_eo.h"
#include "../boundary.h"
#include "../su3.h"
#include "../temporalgauge.h"
#include "../observables.h"
}






int g_numofgpu;

#ifdef GF_8
dev_su3_8 * dev_gf;
dev_su3_8 * h2d_gf;
#else
dev_su3_2v * dev_gf;
dev_su3_2v * h2d_gf;
#endif


dev_spinor* dev_spin1;
dev_spinor* dev_spin2;
dev_spinor* dev_spin3;
dev_spinor* dev_spin4;
dev_spinor* dev_spin5;
dev_spinor* dev_spinin;
dev_spinor* dev_spinout;
dev_spinor * h2d_spin;

//additional spinors for even-odd
dev_spinor* dev_spin_eo1;
dev_spinor* dev_spin_eo2;

int * nn;
int * nn_eo;
int * nn_oe;
int * eoidx_even;
int * eoidx_odd;

int * dev_nn;
int * dev_nn_eo;
int * dev_nn_oe;

int * dev_eoidx_even;
int * dev_eoidx_odd;


size_t output_size;
int* dev_grid;
float * dev_output;
int havedevice = 0;


REAL hostr;
REAL hostkappa;
REAL hostm;
REAL hostmu;


__device__  REAL m;
__device__  REAL mu;
__device__  REAL r=1.0; // this is implicitly assumed to be 1.0 in the host code!!!
__device__  REAL kappa;
__device__ REAL twokappamu;

__device__ dev_complex dev_k0;
__device__ dev_complex dev_k1;
__device__ dev_complex dev_k2;
__device__ dev_complex dev_k3;

__device__ dev_complex dev_mk0;
__device__ dev_complex dev_mk1;
__device__ dev_complex dev_mk2;
__device__ dev_complex dev_mk3;



__constant__ __device__ dev_complex dev_k0c;
__constant__ __device__ dev_complex dev_k1c;
__constant__ __device__ dev_complex dev_k2c;
__constant__ __device__ dev_complex dev_k3c;

__constant__ __device__ dev_complex dev_mk0c;
__constant__ __device__ dev_complex dev_mk1c;
__constant__ __device__ dev_complex dev_mk2c;
__constant__ __device__ dev_complex dev_mk3c;



__device__  int  dev_LX,dev_LY,dev_LZ,dev_T,dev_VOLUME;




////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////




extern "C" {
#include "../Nondegenerate_Matrix.h"
#include "../solver/cg_her_nd.h"
}

#include "MACROS.cuh"



// spinor fields (pointing to device)
dev_spinor * dev_spin1_up;		// auxiliary fields for cg_eo_nd()
dev_spinor * dev_spin1_dn;
dev_spinor * dev_spin2_up;
dev_spinor * dev_spin2_dn;
dev_spinor * dev_spin3_up;
dev_spinor * dev_spin3_dn;
/*
dev_spinor * dev_spin4_up;
dev_spinor * dev_spin4_dn;
dev_spinor * dev_spin5_up;
dev_spinor * dev_spin5_dn;
*/

dev_spinor * dev_spinin_up;		// host/device interaction
dev_spinor * dev_spinin_dn;		// inner/outer interaction
dev_spinor * dev_spinout_up;
dev_spinor * dev_spinout_dn;

dev_spinor * h2d_spin_up;		// for transferring in double precision on host to single precision on device (pointing to host)
dev_spinor * h2d_spin_dn;

dev_spinor * dev_spin_eo1_up;		// auxiliary for  matrix_multiplication32()  called by  dev_cg_eo_nd()
dev_spinor * dev_spin_eo1_dn;
dev_spinor * dev_spin_eo2_up;
dev_spinor * dev_spin_eo2_dn;
dev_spinor * dev_spin_eo3_up;
dev_spinor * dev_spin_eo3_dn;


// physical parameters (on device)
__device__ float mubar, epsbar;


// benchmark
#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
  unsigned long long int device_flops = 0;
#endif

#ifdef CPU_BENCHMARK
  unsigned long long int host_flops = 0;
#endif




////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////




 /* texture for gauge field */
 texture<float4,1, cudaReadModeElementType> gf_tex;
 const textureReference* gf_texRefPtr = NULL;
 cudaChannelFormatDesc gf_channelDesc;
 
 /* texture for spinor field */
 texture<float4,1, cudaReadModeElementType> spin_tex;
 const textureReference* spin_texRefPtr = NULL;
 cudaChannelFormatDesc spin_channelDesc;

 /* texture for spinor field 2*/
 texture<float4,1, cudaReadModeElementType> spin_tex2;
 const textureReference* spin_texRefPtr2 = NULL;
 cudaChannelFormatDesc spin_channelDesc2;


 /* texture for nearest neighbours*/
 texture<int,1, cudaReadModeElementType> nn_tex;
 const textureReference* nn_texRefPtr = NULL;
 cudaChannelFormatDesc nn_channelDesc;



__device__ inline dev_complex dev_cconj (dev_complex c){ /*konjugiert komplexe Zahl*/
 dev_complex erg;
 erg.re = c.re;
 erg.im = -1.0*c.im;
return erg;
}

__device__ inline void dev_ccopy(dev_complex* von, dev_complex* nach){/*kopiert complex von nach complex nach*/
  nach->re = von->re;
  nach->im = von->im;
}

__device__ inline REAL dev_cabssquare (dev_complex c){ /*gibt abs^2 einer komplexen Zahl zurück*/
 return c.re*c.re + c.im*c.im;
}

__device__ inline REAL dev_cabsolute (dev_complex c){/*gibt Betrag einer kompl. zahl zurück*/
 return sqrt(c.re*c.re + c.im*c.im);
}


__device__ inline  dev_complex dev_crealmult(dev_complex c1, REAL real){ /*multipliziert c1 mit reeller zahl re*/
  dev_complex erg;
  erg.re = real*c1.re;
  erg.im = real*c1.im;
return erg;
}

__device__ inline dev_complex dev_cmult (dev_complex c1, dev_complex c2){ /*multiplizier zwei komplexe Zahlen*/
  dev_complex erg;
  erg.re = c1.re * c2.re - c1.im * c2.im;
  erg.im = c1.re * c2.im + c1.im * c2.re;
return erg;
}

__device__ inline dev_complex dev_cadd (dev_complex c1, dev_complex c2){ /*addiert zwei komplexe Zahlen */
  dev_complex erg;
  erg.re = c1.re + c2.re;
  erg.im = c1.im + c2.im;
return erg;
}


__device__ inline dev_complex dev_cdiv(dev_complex c1, dev_complex c2) { /* dividiert c1 durch c2 */
  dev_complex erg;
  REAL oneovernenner = 1.0/(c2.re*c2.re + c2.im*c2.im);
  erg.re = oneovernenner*(c1.re*c2.re + c1.im*c2.im);
  erg.im = oneovernenner*(c1.im*c2.re - c1.re*c2.im);
return erg;
}


__device__ inline dev_complex dev_csub(dev_complex c1, dev_complex c2){
   dev_complex erg;
   erg.re = c1.re - c2.re;
   erg.im = c1.im - c2.im;
return erg;
}


__device__ inline dev_complex dev_initcomplex(REAL re, REAL im){/* gibt komplexe Zahl mit Realt re und Imt im zurück*/
    dev_complex erg;
    erg.re = re;
    erg.im = im;
return (erg);
}





__device__ inline void dev_copy_spinor(dev_spinor *i1, dev_spinor *i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i2+i)).x = (*(i1+i)).x;
    (*(i2+i)).y = (*(i1+i)).y;
    (*(i2+i)).z = (*(i1+i)).z;
    (*(i2+i)).w = (*(i1+i)).w;
  }
}

__device__ inline void dev_zero_spinor(dev_spinor *sin){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(sin+i)).x = 0.0;
    (*(sin+i)).y = 0.0;
    (*(sin+i)).z = 0.0;
    (*(sin+i)).w = 0.0;
  }
}






//out = in + lambda in2
__device__ inline void dev_skalarmult_add_assign_spinor(dev_spinor *in, REAL lambda,dev_spinor * in2, dev_spinor * out){
  int i; 
  #pragma unroll 6
for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + lambda* (*(in2+i)).x;
    (*(out+i)).y = (*(in+i)).y + lambda* (*(in2+i)).y;
    (*(out+i)).z = (*(in+i)).z + lambda* (*(in2+i)).z;
    (*(out+i)).w = (*(in+i)).w + lambda* (*(in2+i)).w;
  }
}




//out = in + lambda in2
__device__ inline void dev_complexmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re - (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + ((*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re - (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + ((*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}




//out = in + (lambda)* in2
__device__ inline void dev_complexcgmult_add_assign_spinor(dev_spinor * in, dev_complex lambda,dev_spinor * in2, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(out+i)).x = (*(in+i)).x + ((*(in2+i)).x*lambda.re + (*(in2+i)).y*lambda.im);
    (*(out+i)).y = (*(in+i)).y + (-(*(in2+i)).x*lambda.im + (*(in2+i)).y*lambda.re);
    (*(out+i)).z = (*(in+i)).z + ((*(in2+i)).z*lambda.re + (*(in2+i)).w*lambda.im);
    (*(out+i)).w = (*(in+i)).w + (-(*(in2+i)).z*lambda.im + (*(in2+i)).w*lambda.re);
  }
}



__device__ void inline dev_skalarmult_spinor(dev_spinor * in, dev_complex lambda, dev_spinor * out){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //out[i] = dev_cmult(in[i],lambda);
    
    (*(out+i)).x = (*(in+i)).x*lambda.re - (*(in+i)).y*lambda.im;
    (*(out+i)).y = (*(in+i)).y*lambda.re + (*(in+i)).x*lambda.im;
    
    (*(out+i)).z = (*(in+i)).z*lambda.re - (*(in+i)).w*lambda.im;
    (*(out+i)).w = (*(in+i)).w*lambda.re + (*(in+i)).z*lambda.im;
  }
}



/*
__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, const dev_complex lambda, dev_spinor * in){


 (*(out)).x = (*(in)).x*lambda.re;
 (*(out)).x -= (*(in)).y*lambda.im;
 
 (*(out)).y = (*(in)).y*lambda.re;
 (*(out)).y += (*(in)).x*lambda.im;
 
 (*(out)).z = (*(in)).z*lambda.re;
 (*(out)).z -= (*(in)).w*lambda.im;

 (*(out)).w = (*(in)).w*lambda.re;
 (*(out)).w += (*(in)).z*lambda.im;

 
 (*(out+1)).x = (*(in+1)).x*lambda.re;
 (*(out+1)).x -= (*(in+1)).y*lambda.im;
 
 (*(out+1)).y = (*(in+1)).y*lambda.re;
 (*(out+1)).y += (*(in+1)).x*lambda.im;
 
 (*(out+1)).z = (*(in+1)).z*lambda.re;
 (*(out+1)).z -= (*(in+1)).w*lambda.im;

 (*(out+1)).w = (*(in+1)).w*lambda.re;
 (*(out+1)).w += (*(in+1)).z*lambda.im;


 (*(out+2)).x = (*(in+2)).x*lambda.re;
 (*(out+2)).x -= (*(in+2)).y*lambda.im;

 (*(out+2)).y = (*(in+2)).y*lambda.re;
 (*(out+2)).y += (*(in+2)).x*lambda.im;
 
 (*(out+2)).z = (*(in+2)).z*lambda.re;
 (*(out+2)).z -= (*(in+2)).w*lambda.im;

 (*(out+2)).w = (*(in+2)).w*lambda.re;
 (*(out+2)).w += (*(in+2)).z*lambda.im;


 (*(out+3)).x = (*(in+3)).y*lambda.im;
 (*(out+3)).x -= (*(in+3)).x*lambda.re;

 (*(out+3)).y = - (*(in+3)).x*lambda.im;
 (*(out+3)).y -= (*(in+3)).y*lambda.re;
 
 (*(out+3)).z = (*(in+3)).w*lambda.im;
 (*(out+3)).z -= (*(in+3)).z*lambda.re;

 (*(out+3)).w = -(*(in+3)).z*lambda.im;
 (*(out+3)).w -= (*(in+3)).w*lambda.re;


 (*(out+4)).x = (*(in+4)).y*lambda.im;
 (*(out+4)).x -= (*(in+4)).x*lambda.re;

 (*(out+4)).y = - (*(in+4)).x*lambda.im;
 (*(out+4)).y -= (*(in+4)).y*lambda.re;
 
 (*(out+4)).z = (*(in+4)).w*lambda.im;
 (*(out+4)).z -= (*(in+4)).z*lambda.re;

 (*(out+4)).w = -(*(in+4)).z*lambda.im;
 (*(out+4)).w -= (*(in+4)).w*lambda.re;


 (*(out+5)).x = (*(in+5)).y*lambda.im;
 (*(out+5)).x -= (*(in+5)).x*lambda.re;

 (*(out+5)).y = - (*(in+5)).x*lambda.im;
 (*(out+5)).y -= (*(in+5)).y*lambda.re;
 
 (*(out+5)).z = (*(in+5)).w*lambda.im;
 (*(out+5)).z -= (*(in+5)).z*lambda.re;

 (*(out+5)).w = -(*(in+5)).z*lambda.im;
 (*(out+5)).w -= (*(in+5)).w*lambda.re;

}
*/


__device__ void inline dev_skalarmult_gamma5_spinor(dev_spinor * out, dev_complex lambda, dev_spinor * in){
int i;
 dev_spinor shelp, tempout;

shelp = *(in);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out)) = tempout;
 
 
shelp = *(in+1);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;
 
 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+1)) = tempout; 
 
 
shelp = *(in+2);
 tempout.x = shelp.x*lambda.re;
 tempout.x -= shelp.y*lambda.im;

 tempout.y = shelp.y*lambda.re;
 tempout.y += shelp.x*lambda.im;
 
 tempout.z = shelp.z*lambda.re;
 tempout.z -= shelp.w*lambda.im;

 tempout.w = shelp.w*lambda.re;
 tempout.w += shelp.z*lambda.im;
(*(out+2)) = tempout; 

 
shelp = *(in+3);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+3)) = tempout;

shelp = *(in+4);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+4)) = tempout;

shelp = *(in+5);
 tempout.x = shelp.y*lambda.im;
 tempout.x -= shelp.x*lambda.re;

 tempout.y = - shelp.x*lambda.im;
 tempout.y -= shelp.y*lambda.re;
 
 tempout.z = shelp.w*lambda.im;
 tempout.z -= shelp.z*lambda.re;

 tempout.w = -shelp.z*lambda.im;
 tempout.w -= shelp.w*lambda.re;
(*(out+5)) = tempout;
}



__device__ void inline dev_realmult_spinor(dev_spinor * in, REAL lambda){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    //in[i] = in[i]*lambda;
    (*(in+i)).x = (*(in+i)).x*lambda;
    (*(in+i)).y = (*(in+i)).y*lambda;
    
    (*(in+i)).z = (*(in+i)).z*lambda;
    (*(in+i)).w = (*(in+i)).w*lambda;
  }
}


__device__ void inline dev_realmult_spinor_assign(dev_spinor* out, REAL lambda, dev_spinor* in){
int i;
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
  //out[i] = in[i]*lambda;
      (*(out+i)).x = (*(in+i)).x*lambda;
	  (*(out+i)).y = (*(in+i)).y*lambda;
      
	  (*(out+i)).z = (*(in+i)).z*lambda;
      (*(out+i)).w = (*(in+i)).w*lambda;
  }
}




__device__ void dev_assign_realmult_add_spinor(dev_spinor* out, REAL lambda, dev_spinor* in1,  dev_spinor* in2){
int i;
REAL help;
//out = lambda*(in1 + in2)
#pragma unroll 6
  for(i=0;i<6;i++){ //color + spin

      help = (*(in1+i)).x*lambda;
      help += (*(in2+i)).x*lambda;
      (*(out+i)).x = help;
      
      help = (*(in1+i)).y*lambda;
      help += (*(in2+i)).y*lambda;
      (*(out+i)).y = help;      

      help = (*(in1+i)).z*lambda;
      help += (*(in2+i)).z*lambda;
      (*(out+i)).z = help;

      help = (*(in1+i)).w*lambda;
      help += (*(in2+i)).w*lambda;
      (*(out+i)).w = help;
  } 
}


__device__ inline void dev_add_spinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x + (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y + (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z + (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w + (*(i2+i)).w;
  }
}



__device__ inline void dev_sub_spinor_assign(dev_spinor * i1, dev_spinor * i2){
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    (*(i1+i)).x = (*(i1+i)).x - (*(i2+i)).x;
    (*(i1+i)).y = (*(i1+i)).y - (*(i2+i)).y;
    (*(i1+i)).z = (*(i1+i)).z - (*(i2+i)).z;
    (*(i1+i)).w = (*(i1+i)).w - (*(i2+i)).w;
  }
}




/*
//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

s1 = tex1Dfetch(spin_tex,6*pos);
s2 = tex1Dfetch(spin_tex,6*pos+1);


//(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
//(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );

(*(out+0)).x =  M[0][0].re*s1.x;
  (*(out+0)).y =  M[0][0].re*s1.y;
(*(out+0)).x -= M[0][0].im*s1.y;
  (*(out+0)).y += M[0][0].im*s1.x;
(*(out+0)).x += M[0][1].re*s1.z;
  (*(out+0)).y += M[0][1].re*s1.w;
(*(out+0)).x -= M[0][1].im*s1.w;
  (*(out+0)).y += M[0][1].im*s1.z;
(*(out+0)).x += M[0][2].re*s2.x;
  (*(out+0)).y += M[0][2].re*s2.y;
(*(out+0)).x -= M[0][2].im*s2.y;
  (*(out+0)).y += M[0][2].im*s2.x;



//(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
//(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+0)).z = M[1][0].re*s1.x;
  (*(out+0)).w = M[1][0].re*s1.y;
(*(out+0)).z -= M[1][0].im*s1.y;
  (*(out+0)).w += M[1][0].im*s1.x;
(*(out+0)).z += M[1][1].re*s1.z;
  (*(out+0)).w += M[1][1].re*s1.w;
(*(out+0)).z -= M[1][1].im*s1.w;
  (*(out+0)).w += M[1][1].im*s1.z;
(*(out+0)).z += M[1][2].re*s2.x;
  (*(out+0)).w += M[1][2].re*s2.y;
(*(out+0)).z -= M[1][2].im*s2.y;
  (*(out+0)).w += M[1][2].im*s2.x;



//(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
//(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );


(*(out+1)).x = M[2][0].re*s1.x;
  (*(out+1)).y = M[2][0].re*s1.y;
(*(out+1)).x -= M[2][0].im*s1.y;
  (*(out+1)).y += M[2][0].im*s1.x;
(*(out+1)).x += M[2][1].re*s1.z;
  (*(out+1)).y += M[2][1].re*s1.w;
(*(out+1)).x -= M[2][1].im*s1.w;
  (*(out+1)).y += M[2][1].im*s1.z;
(*(out+1)).x += M[2][2].re*s2.x;
  (*(out+1)).y += M[2][2].re*s2.y;
(*(out+1)).x -= M[2][2].im*s2.y;
  (*(out+1)).y += M[2][2].im*s2.x;





s1 = tex1Dfetch(spin_tex,6*pos+2);
(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );



s1 = tex1Dfetch(spin_tex,6*pos+3);
s2 = tex1Dfetch(spin_tex,6*pos+4);
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+5);
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}
*/





//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV_spintex(dev_su3 M, int pos, dev_spinor * out){

dev_spinor s1, s2;

s1 = tex1Dfetch(spin_tex,6*pos);
s2 = tex1Dfetch(spin_tex,6*pos+1);


(*(out+0)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+0)).y = ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


//checked by look of eye
(*(out+0)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
// checked 

(*(out+0)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+1)).x = ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+1)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+2);
(*(out+1)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+1)).w =  ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+2)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+2)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+2)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+2)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );



s1 = tex1Dfetch(spin_tex,6*pos+3);
s2 = tex1Dfetch(spin_tex,6*pos+4);
(*(out+3)).x =  ( M[0][0].re*s1.x - M[0][0].im*s1.y ) + ( M[0][1].re*s1.z - M[0][1].im*s1.w ) + ( M[0][2].re*s2.x - M[0][2].im*s2.y );
(*(out+3)).y =   ( M[0][0].re*s1.y + M[0][0].im*s1.x ) + ( M[0][1].re*s1.w + M[0][1].im*s1.z ) + ( M[0][2].re*s2.y + M[0][2].im*s2.x );


(*(out+3)).z =  ( M[1][0].re*s1.x - M[1][0].im*s1.y ) + ( M[1][1].re*s1.z - M[1][1].im*s1.w ) + ( M[1][2].re*s2.x - M[1][2].im*s2.y );
(*(out+3)).w =  ( M[1][0].re*s1.y + M[1][0].im*s1.x ) + ( M[1][1].re*s1.w + M[1][1].im*s1.z ) + ( M[1][2].re*s2.y + M[1][2].im*s2.x );


(*(out+4)).x =  ( M[2][0].re*s1.x - M[2][0].im*s1.y ) + ( M[2][1].re*s1.z - M[2][1].im*s1.w ) + ( M[2][2].re*s2.x - M[2][2].im*s2.y );
(*(out+4)).y =  ( M[2][0].re*s1.y + M[2][0].im*s1.x ) + ( M[2][1].re*s1.w + M[2][1].im*s1.z ) + ( M[2][2].re*s2.y + M[2][2].im*s2.x );



s1 = tex1Dfetch(spin_tex,6*pos+5);
(*(out+4)).z =  ( M[0][0].re*s2.z - M[0][0].im*s2.w ) + ( M[0][1].re*s1.x - M[0][1].im*s1.y ) + ( M[0][2].re*s1.z - M[0][2].im*s1.w );
(*(out+4)).w =   ( M[0][0].re*s2.w + M[0][0].im*s2.z ) + ( M[0][1].re*s1.y + M[0][1].im*s1.x ) + ( M[0][2].re*s1.w + M[0][2].im*s1.z );


(*(out+5)).x = ( M[1][0].re*s2.z - M[1][0].im*s2.w ) + ( M[1][1].re*s1.x - M[1][1].im*s1.y ) + ( M[1][2].re*s1.z - M[1][2].im*s1.w );
(*(out+5)).y =  ( M[1][0].re*s2.w + M[1][0].im*s2.z ) + ( M[1][1].re*s1.y + M[1][1].im*s1.x ) + ( M[1][2].re*s1.w + M[1][2].im*s1.z );


(*(out+5)).z =  ( M[2][0].re*s2.z - M[2][0].im*s2.w ) + ( M[2][1].re*s1.x - M[2][1].im*s1.y ) + ( M[2][2].re*s1.z - M[2][2].im*s1.w );
(*(out+5)).w =  ( M[2][0].re*s2.w + M[2][0].im*s2.z ) + ( M[2][1].re*s1.y + M[2][1].im*s1.x ) + ( M[2][2].re*s1.w + M[2][2].im*s1.z );


}










//multipliziert su3-Matrix mal Spinor im Dirac-Raum
//code in su3_MtV.txt -- generated with codegen
__device__ void dev_su3MtV(dev_su3 M, const dev_spinor * s, dev_spinor * out){

(*(out+0)).x =  ( M[0][0].re*(*(s+0)).x - M[0][0].im*(*(s+0)).y ) + ( M[0][1].re*(*(s+0)).z - M[0][1].im*(*(s+0)).w ) + ( M[0][2].re*(*(s+1)).x - M[0][2].im*(*(s+1)).y );
(*(out+0)).y = ( M[0][0].re*(*(s+0)).y + M[0][0].im*(*(s+0)).x ) + ( M[0][1].re*(*(s+0)).w + M[0][1].im*(*(s+0)).z ) + ( M[0][2].re*(*(s+1)).y + M[0][2].im*(*(s+1)).x );


(*(out+0)).z =  ( M[1][0].re*(*(s+0)).x - M[1][0].im*(*(s+0)).y ) + ( M[1][1].re*(*(s+0)).z - M[1][1].im*(*(s+0)).w ) + ( M[1][2].re*(*(s+1)).x - M[1][2].im*(*(s+1)).y );
(*(out+0)).w =  ( M[1][0].re*(*(s+0)).y + M[1][0].im*(*(s+0)).x ) + ( M[1][1].re*(*(s+0)).w + M[1][1].im*(*(s+0)).z ) + ( M[1][2].re*(*(s+1)).y + M[1][2].im*(*(s+1)).x );


(*(out+1)).x = ( M[2][0].re*(*(s+0)).x - M[2][0].im*(*(s+0)).y ) + ( M[2][1].re*(*(s+0)).z - M[2][1].im*(*(s+0)).w ) + ( M[2][2].re*(*(s+1)).x - M[2][2].im*(*(s+1)).y );
(*(out+1)).y =  ( M[2][0].re*(*(s+0)).y + M[2][0].im*(*(s+0)).x ) + ( M[2][1].re*(*(s+0)).w + M[2][1].im*(*(s+0)).z ) + ( M[2][2].re*(*(s+1)).y + M[2][2].im*(*(s+1)).x );


(*(out+1)).z =  ( M[0][0].re*(*(s+1)).z - M[0][0].im*(*(s+1)).w ) + ( M[0][1].re*(*(s+2)).x - M[0][1].im*(*(s+2)).y ) + ( M[0][2].re*(*(s+2)).z - M[0][2].im*(*(s+2)).w );
(*(out+1)).w =  ( M[0][0].re*(*(s+1)).w + M[0][0].im*(*(s+1)).z ) + ( M[0][1].re*(*(s+2)).y + M[0][1].im*(*(s+2)).x ) + ( M[0][2].re*(*(s+2)).w + M[0][2].im*(*(s+2)).z );


(*(out+2)).x = ( M[1][0].re*(*(s+1)).z - M[1][0].im*(*(s+1)).w ) + ( M[1][1].re*(*(s+2)).x - M[1][1].im*(*(s+2)).y ) + ( M[1][2].re*(*(s+2)).z - M[1][2].im*(*(s+2)).w );
(*(out+2)).y =  ( M[1][0].re*(*(s+1)).w + M[1][0].im*(*(s+1)).z ) + ( M[1][1].re*(*(s+2)).y + M[1][1].im*(*(s+2)).x ) + ( M[1][2].re*(*(s+2)).w + M[1][2].im*(*(s+2)).z );


(*(out+2)).z =  ( M[2][0].re*(*(s+1)).z - M[2][0].im*(*(s+1)).w ) + ( M[2][1].re*(*(s+2)).x - M[2][1].im*(*(s+2)).y ) + ( M[2][2].re*(*(s+2)).z - M[2][2].im*(*(s+2)).w );
(*(out+2)).w =  ( M[2][0].re*(*(s+1)).w + M[2][0].im*(*(s+1)).z ) + ( M[2][1].re*(*(s+2)).y + M[2][1].im*(*(s+2)).x ) + ( M[2][2].re*(*(s+2)).w + M[2][2].im*(*(s+2)).z );


(*(out+3)).x =  ( M[0][0].re*(*(s+3)).x - M[0][0].im*(*(s+3)).y ) + ( M[0][1].re*(*(s+3)).z - M[0][1].im*(*(s+3)).w ) + ( M[0][2].re*(*(s+4)).x - M[0][2].im*(*(s+4)).y );
(*(out+3)).y =   ( M[0][0].re*(*(s+3)).y + M[0][0].im*(*(s+3)).x ) + ( M[0][1].re*(*(s+3)).w + M[0][1].im*(*(s+3)).z ) + ( M[0][2].re*(*(s+4)).y + M[0][2].im*(*(s+4)).x );


(*(out+3)).z =  ( M[1][0].re*(*(s+3)).x - M[1][0].im*(*(s+3)).y ) + ( M[1][1].re*(*(s+3)).z - M[1][1].im*(*(s+3)).w ) + ( M[1][2].re*(*(s+4)).x - M[1][2].im*(*(s+4)).y );
(*(out+3)).w =  ( M[1][0].re*(*(s+3)).y + M[1][0].im*(*(s+3)).x ) + ( M[1][1].re*(*(s+3)).w + M[1][1].im*(*(s+3)).z ) + ( M[1][2].re*(*(s+4)).y + M[1][2].im*(*(s+4)).x );


(*(out+4)).x =  ( M[2][0].re*(*(s+3)).x - M[2][0].im*(*(s+3)).y ) + ( M[2][1].re*(*(s+3)).z - M[2][1].im*(*(s+3)).w ) + ( M[2][2].re*(*(s+4)).x - M[2][2].im*(*(s+4)).y );
(*(out+4)).y =  ( M[2][0].re*(*(s+3)).y + M[2][0].im*(*(s+3)).x ) + ( M[2][1].re*(*(s+3)).w + M[2][1].im*(*(s+3)).z ) + ( M[2][2].re*(*(s+4)).y + M[2][2].im*(*(s+4)).x );


(*(out+4)).z =  ( M[0][0].re*(*(s+4)).z - M[0][0].im*(*(s+4)).w ) + ( M[0][1].re*(*(s+5)).x - M[0][1].im*(*(s+5)).y ) + ( M[0][2].re*(*(s+5)).z - M[0][2].im*(*(s+5)).w );
(*(out+4)).w =   ( M[0][0].re*(*(s+4)).w + M[0][0].im*(*(s+4)).z ) + ( M[0][1].re*(*(s+5)).y + M[0][1].im*(*(s+5)).x ) + ( M[0][2].re*(*(s+5)).w + M[0][2].im*(*(s+5)).z );


(*(out+5)).x = ( M[1][0].re*(*(s+4)).z - M[1][0].im*(*(s+4)).w ) + ( M[1][1].re*(*(s+5)).x - M[1][1].im*(*(s+5)).y ) + ( M[1][2].re*(*(s+5)).z - M[1][2].im*(*(s+5)).w );
(*(out+5)).y =  ( M[1][0].re*(*(s+4)).w + M[1][0].im*(*(s+4)).z ) + ( M[1][1].re*(*(s+5)).y + M[1][1].im*(*(s+5)).x ) + ( M[1][2].re*(*(s+5)).w + M[1][2].im*(*(s+5)).z );


(*(out+5)).z =  ( M[2][0].re*(*(s+4)).z - M[2][0].im*(*(s+4)).w ) + ( M[2][1].re*(*(s+5)).x - M[2][1].im*(*(s+5)).y ) + ( M[2][2].re*(*(s+5)).z - M[2][2].im*(*(s+5)).w );
(*(out+5)).w =  ( M[2][0].re*(*(s+4)).w + M[2][0].im*(*(s+4)).z ) + ( M[2][1].re*(*(s+5)).y + M[2][1].im*(*(s+5)).x ) + ( M[2][2].re*(*(s+5)).w + M[2][2].im*(*(s+5)).z );
}





//multipliziert gedaggerte su3-Matrix mal Spinor im Dirac-Raum  -- generated with codegen
__device__ void dev_su3MdaggertV(dev_su3 M, dev_spinor * s, dev_spinor * out){
  dev_complex help1;
help1.re = M[0][0].re*(*(s+0)).x + M[0][0].im*(*(s+0)).y + M[1][0].re*(*(s+0)).z + M[1][0].im*(*(s+0)).w + M[2][0].re*(*(s+1)).x + M[2][0].im*(*(s+1)).y;
(*(out+0)).x = help1.re;
help1.im = M[0][0].re*(*(s+0)).y - M[0][0].im*(*(s+0)).x + M[1][0].re*(*(s+0)).w - M[1][0].im*(*(s+0)).z + M[2][0].re*(*(s+1)).y - M[2][0].im*(*(s+1)).x;
(*(out+0)).y = help1.im;

help1.re = M[0][1].re*(*(s+0)).x + M[0][1].im*(*(s+0)).y + M[1][1].re*(*(s+0)).z + M[1][1].im*(*(s+0)).w + M[2][1].re*(*(s+1)).x + M[2][1].im*(*(s+1)).y;
(*(out+0)).z = help1.re;
help1.im = M[0][1].re*(*(s+0)).y - M[0][1].im*(*(s+0)).x + M[1][1].re*(*(s+0)).w - M[1][1].im*(*(s+0)).z + M[2][1].re*(*(s+1)).y - M[2][1].im*(*(s+1)).x;
(*(out+0)).w = help1.im;

help1.re = M[0][2].re*(*(s+0)).x + M[0][2].im*(*(s+0)).y + M[1][2].re*(*(s+0)).z + M[1][2].im*(*(s+0)).w + M[2][2].re*(*(s+1)).x + M[2][2].im*(*(s+1)).y;
(*(out+1)).x = help1.re;
help1.im = M[0][2].re*(*(s+0)).y - M[0][2].im*(*(s+0)).x + M[1][2].re*(*(s+0)).w - M[1][2].im*(*(s+0)).z + M[2][2].re*(*(s+1)).y - M[2][2].im*(*(s+1)).x;
(*(out+1)).y = help1.im;

help1.re = M[0][0].re*(*(s+1)).z + M[0][0].im*(*(s+1)).w + M[1][0].re*(*(s+2)).x + M[1][0].im*(*(s+2)).y + M[2][0].re*(*(s+2)).z + M[2][0].im*(*(s+2)).w;
(*(out+1)).z = help1.re;
help1.im = M[0][0].re*(*(s+1)).w - M[0][0].im*(*(s+1)).z + M[1][0].re*(*(s+2)).y - M[1][0].im*(*(s+2)).x + M[2][0].re*(*(s+2)).w - M[2][0].im*(*(s+2)).z;
(*(out+1)).w = help1.im;

help1.re = M[0][1].re*(*(s+1)).z + M[0][1].im*(*(s+1)).w + M[1][1].re*(*(s+2)).x + M[1][1].im*(*(s+2)).y + M[2][1].re*(*(s+2)).z + M[2][1].im*(*(s+2)).w;
(*(out+2)).x = help1.re;
help1.im = M[0][1].re*(*(s+1)).w - M[0][1].im*(*(s+1)).z + M[1][1].re*(*(s+2)).y - M[1][1].im*(*(s+2)).x + M[2][1].re*(*(s+2)).w - M[2][1].im*(*(s+2)).z;
(*(out+2)).y = help1.im;

help1.re = M[0][2].re*(*(s+1)).z + M[0][2].im*(*(s+1)).w + M[1][2].re*(*(s+2)).x + M[1][2].im*(*(s+2)).y + M[2][2].re*(*(s+2)).z + M[2][2].im*(*(s+2)).w;
(*(out+2)).z = help1.re;
help1.im = M[0][2].re*(*(s+1)).w - M[0][2].im*(*(s+1)).z + M[1][2].re*(*(s+2)).y - M[1][2].im*(*(s+2)).x + M[2][2].re*(*(s+2)).w - M[2][2].im*(*(s+2)).z;
(*(out+2)).w = help1.im;

help1.re = M[0][0].re*(*(s+3)).x + M[0][0].im*(*(s+3)).y + M[1][0].re*(*(s+3)).z + M[1][0].im*(*(s+3)).w + M[2][0].re*(*(s+4)).x + M[2][0].im*(*(s+4)).y;
(*(out+3)).x = help1.re;
help1.im = M[0][0].re*(*(s+3)).y - M[0][0].im*(*(s+3)).x + M[1][0].re*(*(s+3)).w - M[1][0].im*(*(s+3)).z + M[2][0].re*(*(s+4)).y - M[2][0].im*(*(s+4)).x;
(*(out+3)).y = help1.im;

help1.re = M[0][1].re*(*(s+3)).x + M[0][1].im*(*(s+3)).y + M[1][1].re*(*(s+3)).z + M[1][1].im*(*(s+3)).w + M[2][1].re*(*(s+4)).x + M[2][1].im*(*(s+4)).y;
(*(out+3)).z = help1.re;
help1.im = M[0][1].re*(*(s+3)).y - M[0][1].im*(*(s+3)).x + M[1][1].re*(*(s+3)).w - M[1][1].im*(*(s+3)).z + M[2][1].re*(*(s+4)).y - M[2][1].im*(*(s+4)).x;
(*(out+3)).w = help1.im;

help1.re = M[0][2].re*(*(s+3)).x + M[0][2].im*(*(s+3)).y + M[1][2].re*(*(s+3)).z + M[1][2].im*(*(s+3)).w + M[2][2].re*(*(s+4)).x + M[2][2].im*(*(s+4)).y;
(*(out+4)).x = help1.re;
help1.im = M[0][2].re*(*(s+3)).y - M[0][2].im*(*(s+3)).x + M[1][2].re*(*(s+3)).w - M[1][2].im*(*(s+3)).z + M[2][2].re*(*(s+4)).y - M[2][2].im*(*(s+4)).x;
(*(out+4)).y = help1.im;

help1.re = M[0][0].re*(*(s+4)).z + M[0][0].im*(*(s+4)).w + M[1][0].re*(*(s+5)).x + M[1][0].im*(*(s+5)).y + M[2][0].re*(*(s+5)).z + M[2][0].im*(*(s+5)).w;
(*(out+4)).z = help1.re;
help1.im = M[0][0].re*(*(s+4)).w - M[0][0].im*(*(s+4)).z + M[1][0].re*(*(s+5)).y - M[1][0].im*(*(s+5)).x + M[2][0].re*(*(s+5)).w - M[2][0].im*(*(s+5)).z;
(*(out+4)).w = help1.im;

help1.re = M[0][1].re*(*(s+4)).z + M[0][1].im*(*(s+4)).w + M[1][1].re*(*(s+5)).x + M[1][1].im*(*(s+5)).y + M[2][1].re*(*(s+5)).z + M[2][1].im*(*(s+5)).w;
(*(out+5)).x = help1.re;
help1.im = M[0][1].re*(*(s+4)).w - M[0][1].im*(*(s+4)).z + M[1][1].re*(*(s+5)).y - M[1][1].im*(*(s+5)).x + M[2][1].re*(*(s+5)).w - M[2][1].im*(*(s+5)).z;
(*(out+5)).y = help1.im;

help1.re = M[0][2].re*(*(s+4)).z + M[0][2].im*(*(s+4)).w + M[1][2].re*(*(s+5)).x + M[1][2].im*(*(s+5)).y + M[2][2].re*(*(s+5)).z + M[2][2].im*(*(s+5)).w;
(*(out+5)).z = help1.re;
help1.im = M[0][2].re*(*(s+4)).w - M[0][2].im*(*(s+4)).z + M[1][2].re*(*(s+5)).y - M[1][2].im*(*(s+5)).x + M[2][2].re*(*(s+5)).w - M[2][2].im*(*(s+5)).z;
(*(out+5)).w = help1.im;
}




// Gamma t
__device__ void dev_Gamma0(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;
}



//Gamma z
__device__ void dev_Gamma3(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;

}



//Gamma y
__device__ void dev_Gamma2(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;
}



//Gamma x
__device__ void dev_Gamma1(dev_spinor * in){
  REAL tempre,tempim;
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
  
}



__device__ void dev_Gamma5(dev_spinor * in){
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
}


__device__ void dev_Gamma5_assign(dev_spinor* out, dev_spinor* in){
  (*(out)).x = (*(in)).x;
  (*(out)).y = (*(in)).y;
  (*(out)).z = (*(in)).z;
  (*(out)).w = (*(in)).w;
  (*(out+1)).x = (*(in+1)).x;
  (*(out+1)).y = (*(in+1)).y;

  (*(out+1)).z = (*(in+1)).z;
  (*(out+1)).w = (*(in+1)).w;
  (*(out+2)).x = (*(in+2)).x;
  (*(out+2)).y = (*(in+2)).y;
  (*(out+2)).z = (*(in+2)).z;
  (*(out+2)).w = (*(in+2)).w;

  (*(out+3)).x = -1.0*(*(in+3)).x;
  (*(out+3)).y = -1.0*(*(in+3)).y;
  (*(out+3)).z = -1.0*(*(in+3)).z;
  (*(out+3)).w = -1.0*(*(in+3)).w;
  (*(out+4)).x = -1.0*(*(in+4)).x;
  (*(out+4)).y = -1.0*(*(in+4)).y;

  (*(out+4)).z = -1.0*(*(in+4)).z;
  (*(out+4)).w = -1.0*(*(in+4)).w;
  (*(out+5)).x = -1.0*(*(in+5)).x;
  (*(out+5)).y = -1.0*(*(in+5)).y;
  (*(out+5)).z = -1.0*(*(in+5)).z;
  (*(out+5)).w = -1.0*(*(in+5)).w;

}




// older version, all in one function
__device__ void dev_GammatV(int mu, dev_spinor * in){//multipliziert Gamma(mu)*V effizientes ausnutzen der Nullen 
 REAL tempre,tempim;
 /* ORDER: t, z, y, x*/
 switch (mu){
 
 case 0:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+3)).x;
     (*(in+0)).y = -1.0*(*(in+3)).y;
     (*(in+3)).x = -1.0*tempre;
     (*(in+3)).y = -1.0*tempim;     
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+3)).z;
     (*(in+0)).w = -1.0*(*(in+3)).w;
     (*(in+3)).z = -1.0*tempre;
     (*(in+3)).w = -1.0*tempim; 
 
 
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+4)).x;
     (*(in+1)).y = -1.0*(*(in+4)).y;
     (*(in+4)).x = -1.0*tempre;
     (*(in+4)).y = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).z;
     (*(in+1)).w = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).x;
     (*(in+2)).y = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
   
   
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).z;
     (*(in+2)).w = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;

 break;
 
 case 1:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+3)).y;
     (*(in+0)).y = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+3)).w;
     (*(in+0)).w = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;    
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+4)).y;
     (*(in+1)).y = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;     
     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = -1.0*(*(in+4)).w;
     (*(in+1)).w = (*(in+4)).z;
     (*(in+4)).z  = tempim;
     (*(in+4)).w  = -1.0*tempre;     
     
     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = -1.0*(*(in+5)).y;
     (*(in+2)).y = (*(in+5)).x;
     (*(in+5)).x = tempim;
     (*(in+5)).y = -1.0*tempre;    
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = -1.0*(*(in+5)).w;
     (*(in+2)).w = (*(in+5)).z;
     (*(in+5)).z = tempim;
     (*(in+5)).w = -1.0*tempre;


 break;
 
 case 2:
     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = -1.0*(*(in+4)).z;
     (*(in+0)).y = -1.0*(*(in+4)).w;
     (*(in+4)).z = -1.0*tempre;
     (*(in+4)).w = -1.0*tempim;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = -1.0*(*(in+5)).x;
     (*(in+0)).w = -1.0*(*(in+5)).y;
     (*(in+5)).x = -1.0*tempre;
     (*(in+5)).y = -1.0*tempim;     
     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = -1.0*(*(in+5)).z;
     (*(in+1)).y = -1.0*(*(in+5)).w;
     (*(in+5)).z = -1.0*tempre;
     (*(in+5)).w = -1.0*tempim;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).x;
     (*(in+1)).w = (*(in+3)).y;
     (*(in+3)).x = tempre;
     (*(in+3)).y = tempim;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).z;
     (*(in+2)).y = (*(in+3)).w;
     (*(in+3)).z = tempre;
     (*(in+3)).w = tempim;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).x;
     (*(in+2)).w = (*(in+4)).y;
     (*(in+4)).x = tempre;
     (*(in+4)).y = tempim;

 break; 
 
 case 3:


     tempre = (*(in+0)).x;
     tempim = (*(in+0)).y;
     (*(in+0)).x = (*(in+4)).w;
     (*(in+0)).y = -1.0*(*(in+4)).z;
     (*(in+4)).z  = -1.0*tempim;
     (*(in+4)).w  = tempre;    
     
     tempre = (*(in+0)).z;
     tempim = (*(in+0)).w;
     (*(in+0)).z = (*(in+5)).y;
     (*(in+0)).w = -1.0*(*(in+5)).x;
     (*(in+5)).x = -1.0*tempim;
     (*(in+5)).y = tempre;     
     
     tempre = (*(in+1)).x;
     tempim = (*(in+1)).y;
     (*(in+1)).x = (*(in+5)).w;
     (*(in+1)).y = -1.0*(*(in+5)).z;
     (*(in+5)).z = -1.0*tempim;
     (*(in+5)).w = tempre;     
     
     tempre = (*(in+1)).z;
     tempim = (*(in+1)).w;
     (*(in+1)).z = (*(in+3)).y;
     (*(in+1)).w = -1.0*(*(in+3)).x;
     (*(in+3)).x = -1.0*tempim;
     (*(in+3)).y = tempre;     
     
     tempre = (*(in+2)).x;
     tempim = (*(in+2)).y;
     (*(in+2)).x = (*(in+3)).w;
     (*(in+2)).y = -1.0*(*(in+3)).z;
     (*(in+3)).z = -1.0*tempim;
     (*(in+3)).w = tempre;     
     
     
     tempre = (*(in+2)).z;
     tempim = (*(in+2)).w;
     (*(in+2)).z = (*(in+4)).y;
     (*(in+2)).w = -1.0*(*(in+4)).x;
     (*(in+4)).x = -1.0*tempim;
     (*(in+4)).y = tempre;
     
     
 break;
 
 
 case 4:
  
          (*(in+3)).x = -1.0*(*(in+3)).x;
          (*(in+3)).y = -1.0*(*(in+3)).y;
          (*(in+3)).z = -1.0*(*(in+3)).z;
          (*(in+3)).w = -1.0*(*(in+3)).w;
          (*(in+4)).x = -1.0*(*(in+4)).x;
          (*(in+4)).y = -1.0*(*(in+4)).y; 

          (*(in+4)).z = -1.0*(*(in+4)).z;
          (*(in+4)).w = -1.0*(*(in+4)).w;
          (*(in+5)).x = -1.0*(*(in+5)).x;
          (*(in+5)).y = -1.0*(*(in+5)).y;
          (*(in+5)).z = -1.0*(*(in+5)).z;
          (*(in+5)).w = -1.0*(*(in+5)).w;  
 break;
 }
}




// reconstruction of the link fields from two rows of the su3 matrix
// numbers are fetched from texture cache
__device__ void dev_reconstructgf_2vtexref (const dev_su3_2v* field, int pos, dev_su3* gf){
  dev_complex help1;
  dev_complex help2;
  float4 gfin;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = field[3*pos];
  #endif
  //first row
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = gfin.y;
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = field[3*pos + 1];
  #endif
  (*gf)[0][2].re = gfin.x;
  (*gf)[0][2].im = gfin.y;
  //second row
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = field[3*pos + 2];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = gfin.y;
  (*gf)[1][2].re = gfin.z;
  (*gf)[1][2].im = gfin.w;
  
  //third row from cconj(cross product of first and second row)

  (*gf)[2][0].re = (*gf)[0][1].re * (*gf)[1][2].re;
  (*gf)[2][0].re -= (*gf)[0][1].im * (*gf)[1][2].im;
  (*gf)[2][0].re -= (*gf)[0][2].re * (*gf)[1][1].re;
  (*gf)[2][0].re += (*gf)[0][2].im * (*gf)[1][1].im;
  
  (*gf)[2][0].im = -(*gf)[0][1].re * (*gf)[1][2].im;
  (*gf)[2][0].im -= (*gf)[0][1].im * (*gf)[1][2].re;
  (*gf)[2][0].im += (*gf)[0][2].re * (*gf)[1][1].im;
  (*gf)[2][0].im += (*gf)[0][2].im * (*gf)[1][1].re;
  


  (*gf)[2][1].re = (*gf)[0][2].re * (*gf)[1][0].re;
  (*gf)[2][1].re -= (*gf)[0][2].im * (*gf)[1][0].im;
  (*gf)[2][1].re -= (*gf)[0][0].re * (*gf)[1][2].re;
  (*gf)[2][1].re += (*gf)[0][0].im * (*gf)[1][2].im;
  
  (*gf)[2][1].im = -(*gf)[0][2].re * (*gf)[1][0].im;
  (*gf)[2][1].im -= (*gf)[0][2].im * (*gf)[1][0].re;
  (*gf)[2][1].im += (*gf)[0][0].re * (*gf)[1][2].im;
  (*gf)[2][1].im += (*gf)[0][0].im * (*gf)[1][2].re;


  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re;
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im;
  (*gf)[2][2].re -= (*gf)[0][1].re * (*gf)[1][0].re;
  (*gf)[2][2].re += (*gf)[0][1].im * (*gf)[1][0].im;
  
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re;
  (*gf)[2][2].im += (*gf)[0][1].re * (*gf)[1][0].im;
  (*gf)[2][2].im += (*gf)[0][1].im * (*gf)[1][0].re;
  

  return;
}




// su3 - dagger reconstruction from two rows  
__device__ void dev_reconstructgf_2vtexref_dagger (const dev_su3_2v* field, int pos, dev_su3* gf){
  //dev_complex help1;
  //dev_complex help2;
  float4 gfin;
  
  
  //first column (minus in im for complex conj.)
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos);
  #else
    gfin = field[3*pos];
  #endif
  (*gf)[0][0].re = gfin.x;
  (*gf)[0][0].im = -gfin.y;
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+1);
  #else
    gfin = field[3*pos +1];
  #endif
  (*gf)[2][0].re = gfin.x;
  (*gf)[2][0].im = -gfin.y;
  
  //second  column (minus in im for complex conj.)
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,3*pos+2);
  #else
    gfin = field[3*pos +2];
  #endif
  (*gf)[1][1].re = gfin.x;
  (*gf)[1][1].im = -gfin.y;
  (*gf)[2][1].re = gfin.z;
  (*gf)[2][1].im = -gfin.w;
  


  (*gf)[0][2].re = (*gf)[1][0].re * (*gf)[2][1].re;
  (*gf)[0][2].re -= (*gf)[1][0].im * (*gf)[2][1].im;
  (*gf)[0][2].re -= (*gf)[2][0].re * (*gf)[1][1].re;
  (*gf)[0][2].re += (*gf)[2][0].im * (*gf)[1][1].im;
  
  (*gf)[0][2].im = -(*gf)[1][0].re* (*gf)[2][1].im;
  (*gf)[0][2].im -= (*gf)[1][0].im* (*gf)[2][1].re;
  (*gf)[0][2].im += (*gf)[2][0].re*(*gf)[1][1].im;
  (*gf)[0][2].im += (*gf)[2][0].im*(*gf)[1][1].re;
  

  (*gf)[1][2].re = (*gf)[2][0].re*(*gf)[0][1].re;
  (*gf)[1][2].re -= (*gf)[2][0].im*(*gf)[0][1].im;
  (*gf)[1][2].re -= (*gf)[0][0].re*(*gf)[2][1].re;
  (*gf)[1][2].re += (*gf)[0][0].im*(*gf)[2][1].im;
  
  (*gf)[1][2].im = -(*gf)[2][0].re * (*gf)[0][1].im;
  (*gf)[1][2].im -= (*gf)[2][0].im * (*gf)[0][1].re;
  (*gf)[1][2].im += (*gf)[0][0].re * (*gf)[2][1].im;
  (*gf)[1][2].im += (*gf)[0][0].im * (*gf)[2][1].re;
  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[1][1].re;
  (*gf)[2][2].re -= (*gf)[0][0].im * (*gf)[1][1].im;
  (*gf)[2][2].re -= (*gf)[1][0].re * (*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[1][0].im * (*gf)[0][1].im;
  
  (*gf)[2][2].im = -(*gf)[0][0].re * (*gf)[1][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[1][1].re;
  (*gf)[2][2].im += (*gf)[1][0].re * (*gf)[0][1].im;
  (*gf)[2][2].im += (*gf)[1][0].im * (*gf)[0][1].re;
  
}






// reconstruction of the gf using 8 real parameters as 
// described in the appendix of hep-lat 0911.3191 (M.Clark et al.)
// optimized once
__device__ void dev_reconstructgf_8texref (const dev_su3_2v * field, int pos, dev_su3* gf){

  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = field[2*pos];
  #endif
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re); //reciprocal sqrt

  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = field[2*pos + 1];
  #endif
  
  // reconstruct a1 use sqrt instead of sin
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrtf(help);
  }
  else{
    p1.re = 0.0f;
  }
  
  sincos(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = (*gf)[0][0].im * p1.re;
  
  
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1 use sqrt instead of sin
  help =1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im;
  if(help > 0.0f){
    p1.re = sqrtf(help);
  }   
  else{      
    p1.re = 0.0f;  
  }
  sincos(gfin.y, &(*gf)[2][0].im, &(*gf)[2][0].re);
  (*gf)[2][0].re = (*gf)[2][0].re * p1.re; 
  (*gf)[2][0].im = (*gf)[2][0].im * p1.re;
   
  
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  
  (*gf)[1][1].re = p1.re*(*gf)[0][2].re;
  (*gf)[1][1].re += p1.im*(*gf)[0][2].im;
  (*gf)[1][1].im = p1.im*(*gf)[0][2].re;
  (*gf)[1][1].im -= p1.re*(*gf)[0][2].im;
  
  (*gf)[2][2].re = (*gf)[0][0].re * (*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[0][0].im * (*gf)[0][1].im;
  
  (*gf)[2][2].im = (*gf)[0][0].re * (*gf)[0][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im * (*gf)[0][1].re;
  (*gf)[2][2] = dev_cmult(p2, (*gf)[2][2]);
  
  (*gf)[1][1].re = -one_over_N*( (*gf)[1][1].re + (*gf)[2][2].re);
  (*gf)[1][1].im = -one_over_N*((*gf)[1][1].im + (*gf)[2][2].im);
  
  
  
  
  
  // calculate b3
  (*gf)[1][2].re = p1.re*(*gf)[0][1].re;
  (*gf)[1][2].re += p1.im*(*gf)[0][1].im;
  (*gf)[1][2].im = p1.im*(*gf)[0][1].re;
  (*gf)[1][2].im -= p1.re*(*gf)[0][1].im;
  
  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][2].re;
  (*gf)[2][2].re += (*gf)[0][0].im*(*gf)[0][2].im;
  (*gf)[2][2].im = (*gf)[0][0].re*(*gf)[0][2].im;
  (*gf)[2][2].im -= (*gf)[0][0].im*(*gf)[0][2].re;
  (*gf)[2][2] = dev_cmult(p2,(*gf)[2][2]);
  
  (*gf)[1][2].re = one_over_N*( (*gf)[1][2].re - (*gf)[2][2].re);
  (*gf)[1][2].im = one_over_N*( (*gf)[1][2].im - (*gf)[2][2].im);
  
  
  // calculate c2
  (*gf)[2][1].re = p2.re*(*gf)[0][2].re;
  (*gf)[2][1].re -= p2.im*(*gf)[0][2].im;
  (*gf)[2][1].im = -p2.re*(*gf)[0][2].im;
  (*gf)[2][1].im -= p2.im*(*gf)[0][2].re;
  
  

  (*gf)[2][2].re = (*gf)[0][0].re*(*gf)[0][1].re;
  (*gf)[2][2].re += (*gf)[0][0].im*(*gf)[0][1].im;
  (*gf)[2][2].im = (*gf)[0][0].re* (*gf)[0][1].im;
  (*gf)[2][2].im -= (*gf)[0][0].im* (*gf)[0][1].re;
  help = (*gf)[2][2].re;
  (*gf)[2][2].re = p1.re*(*gf)[2][2].re;
  (*gf)[2][2].re += p1.im*(*gf)[2][2].im;
  (*gf)[2][2].im = p1.re*(*gf)[2][2].im - p1.im*help;
  
  
  (*gf)[2][1].re = one_over_N*((*gf)[2][1].re - (*gf)[2][2].re);
  (*gf)[2][1].im = one_over_N*((*gf)[2][1].im - (*gf)[2][2].im);
  
  // now we have to use p2 and p1 as a help variable, as this is not 
  // needed any more after the first
  // step
  // calculate c3
  (*gf)[2][2].re = p2.re * (*gf)[0][1].re;
  (*gf)[2][2].re -= p2.im * (*gf)[0][1].im;
  (*gf)[2][2].im = - p2.im*(*gf)[0][1].re;
  (*gf)[2][2].im -= p2.re*(*gf)[0][1].im;
  
  p2.re = (*gf)[0][0].re * (*gf)[0][2].re;
  p2.re += (*gf)[0][0].im * (*gf)[0][2].im;
  p2.im = (*gf)[0][0].re * (*gf)[0][2].im;
  p2.im -= (*gf)[0][0].im * (*gf)[0][2].re;
  p2 = dev_cmult(  dev_cconj(p1) , p2);
  
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_crealmult((*gf)[2][2], -one_over_N);
                      
}








__device__ void dev_reconstructgf_8texref_dagger (const dev_su3_2v* field,int pos, dev_su3* gf){


  float4 gfin;
  REAL one_over_N, help;
  dev_complex p1,p2;
  
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos);
  #else
    gfin = field[2*pos];
  #endif
  // read a2 a3
  (*gf)[1][0].re = gfin.x;
  (*gf)[1][0].im = -gfin.y;
  (*gf)[2][0].re = gfin.z;
  (*gf)[2][0].im = -gfin.w;  
 
  p1.re = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  one_over_N = rsqrtf(p1.re);  // reciprocal sqrt

  
  // read theta_a1, theta_c1, b1
  #ifdef USETEXTURE
    gfin = tex1Dfetch(gf_tex,2*pos + 1);
  #else
    gfin = field[2*pos + 1];
  #endif
  
  // reconstruct a1
  help = 1.0f - p1.re;
  if(help > 0.0f){
     p1.re = sqrtf(help);   
  }
  else{
    p1.re = 0.0f;
  }
  //(*gf)[0][0].re = p1.re*cosf(gfin.x);
  //(*gf)[0][0].im = -p1.re*sinf(gfin.x);
  
  sincos(gfin.x, &(*gf)[0][0].im, &(*gf)[0][0].re);
  (*gf)[0][0].re = (*gf)[0][0].re * p1.re;
  (*gf)[0][0].im = -(*gf)[0][0].im * p1.re;
    
  
  
  
  // assign b1
  (*gf)[0][1].re = gfin.z;
  (*gf)[0][1].im = -gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[0][1].re;
  p2.im = -one_over_N*(*gf)[0][1].im;  


  // reconstruct c1
  help = 1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[0][1].re * (*gf)[0][1].re - (*gf)[0][1].im * (*gf)[0][1].im;
  if(help > 0.0f){
    p1.re = sqrtf(help);
  }
  else{
    p1.re = 0.0f;
  }
  //(*gf)[0][2].re = p1.re*cosf(gfin.y);
  //(*gf)[0][2].im = -p1.re*sinf(gfin.y);
  
  sincos(gfin.y, &(*gf)[0][2].im, &(*gf)[0][2].re);
  (*gf)[0][2].re = (*gf)[0][2].re * p1.re;
  (*gf)[0][2].im = -(*gf)[0][2].im * p1.re;
     
  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[0][2].re;
  p1.im = one_over_N*(*gf)[0][2].im;
  
  //use the last reconstructed gf component gf[2][2] (c3) as a help variable for b2,b3 and c2
  //this is in order to save registers and to prevent extra loading and storing from global mem
  // calculate b2
  (*gf)[1][1] = dev_cmult(p1,   (*gf)[2][0]    );
  (*gf)[2][2] = dev_cmult(p2, dev_cmult( (*gf)[0][0] , dev_cconj((*gf)[1][0] ))  );
  (*gf)[1][1] = dev_cadd((*gf)[1][1], (*gf)[2][2]);
  (*gf)[1][1] = dev_cconj(dev_crealmult((*gf)[1][1], -one_over_N));
  
  // calculate b3
  (*gf)[2][1] = dev_cmult(p1,   (*gf)[1][0]    );
  (*gf)[2][2] = dev_cmult(p2, dev_cmult( (*gf)[0][0] , dev_cconj((*gf)[2][0] ))  );
  (*gf)[2][1] = dev_csub((*gf)[2][1], (*gf)[2][2]);
  (*gf)[2][1] = dev_cconj(dev_crealmult((*gf)[2][1], one_over_N));
  
  // calculate c2
  (*gf)[1][2] = dev_cmult(  dev_cconj(p2) ,  (*gf)[2][0]    );
  (*gf)[2][2] = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(   (*gf)[0][0]  , dev_cconj( (*gf)[1][0]) )
                     );
  (*gf)[1][2] = dev_csub((*gf)[1][2], (*gf)[2][2]);
  (*gf)[1][2] = dev_cconj(dev_crealmult((*gf)[1][2], one_over_N));
  
  // use p2 as help variable after the first step
  // calculate c3
  (*gf)[2][2] = dev_cmult(  dev_cconj(p2) ,   (*gf)[1][0]    );
  p2 = dev_cmult(  dev_cconj(p1) , 
                       dev_cmult(   (*gf)[0][0]  , dev_cconj((*gf)[2][0] ) )
                     );
  (*gf)[2][2] = dev_cadd((*gf)[2][2], p2);
  (*gf)[2][2] = dev_cconj(dev_crealmult((*gf)[2][2], -one_over_N));

}





__global__ void dev_gamma5(dev_spinor * sin, dev_spinor * sout){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          sout[6*pos+0].x = sin[6*pos+0].x;
          sout[6*pos+0].y = sin[6*pos+0].y;
          sout[6*pos+0].z = sin[6*pos+0].z;
          sout[6*pos+0].w = sin[6*pos+0].w;
          sout[6*pos+1].x = sin[6*pos+1].x;
          sout[6*pos+1].y = sin[6*pos+1].y;
          
          sout[6*pos+1].z = sin[6*pos+1].z;
          sout[6*pos+1].w = sin[6*pos+1].w;
          sout[6*pos+2].x = sin[6*pos+2].x;
          sout[6*pos+2].y = sin[6*pos+2].y;
          sout[6*pos+2].z = sin[6*pos+2].z;
          sout[6*pos+2].w = sin[6*pos+2].w;   
          
          sout[6*pos+3].x = -1.0*sin[6*pos+3].x;
          sout[6*pos+3].y = -1.0*sin[6*pos+3].y;
          sout[6*pos+3].z = -1.0*sin[6*pos+3].z;
          sout[6*pos+3].w = -1.0*sin[6*pos+3].w;
          sout[6*pos+4].x = -1.0*sin[6*pos+4].x;
          sout[6*pos+4].y = -1.0*sin[6*pos+4].y; 

          sout[6*pos+4].z = -1.0*sin[6*pos+4].z;
          sout[6*pos+4].w = -1.0*sin[6*pos+4].w;
          sout[6*pos+5].x = -1.0*sin[6*pos+5].x;
          sout[6*pos+5].y = -1.0*sin[6*pos+5].y;
          sout[6*pos+5].z = -1.0*sin[6*pos+5].z;
          sout[6*pos+5].w = -1.0*sin[6*pos+5].w;                 
  } 
}




__global__ void dev_swapmu(){
  if(blockIdx.x == 0 && threadIdx.x == 0){
    mu = - mu;
  }
}






// computes sout = 1/(1 +- mutilde gamma5) sin = (1 -+ i mutilde gamma5)/(1+mutilde^2) sin
// mutilde = 2 kappa mu
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_inv(dev_spinor* sin, dev_spinor* sout, const REAL sign){
   
   dev_spinor slocal[6];
   //need the inverse sign in the numerator because of inverse
   dev_complex pm_imu = dev_initcomplex(0.0,-1.0*sign*twokappamu);
   
   REAL one_plus_musquare_inv = 1.0/(1.0 + twokappamu*twokappamu);
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x;  
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );
	 dev_add_spinor_assign(&(slocal[0]), &(sin[6*pos]));
     //dev_realmult_spinor(&(slocal[0]), one_plus_musquare_inv);
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos])); 
     dev_realmult_spinor_assign(&(sout[6*pos]), one_plus_musquare_inv, &(slocal[0]) );
   }
}





// sout = gamma_5*((1\pm i\mutilde \gamma_5)*sin1 - sin2)
// uses shared local memory for manipulation
__global__ void dev_mul_one_pm_imu_sub_mul_gamma5(dev_spinor* sin1, dev_spinor* sin2, dev_spinor* sout, const REAL sign){
   dev_spinor slocal[6];
   dev_complex pm_imu = dev_initcomplex(0.0, sign*twokappamu); // i mutilde
   int pos;
   pos= threadIdx.x + blockDim.x*blockIdx.x; 
   int ix = threadIdx.x;
   if(pos < dev_VOLUME){
     //dev_skalarmult_spinor(&(sin1[6*pos]), pm_imu, &(slocal[0]));
     //dev_Gamma5(&(slocal[0]));
     dev_skalarmult_gamma5_spinor(&(slocal[0]),pm_imu,&(sin1[6*pos]));
	 dev_add_spinor_assign(&(slocal[0]), &(sin1[6*pos]));
     dev_sub_spinor_assign(&(slocal[0]), &(sin2[6*pos]));
     //dev_Gamma5(&(slocal[0]));
     //dev_copy_spinor(&(slocal[0]), &(sout[6*pos]));
     dev_Gamma5_assign(&(sout[6*pos]), &(slocal[0]));
   }   
}






//-kappa(r - gamma_mu)
__device__ void dev_kappaP1_plus(dev_spinor * out, dev_spinor * in, REAL kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+4)).w);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+4)).z);
     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+5)).y);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+5)).x);    

     (*(out+1)).x -= kappa*((*(in+1)).x - (*(in+5)).w);
     (*(out+1)).y -= kappa*((*(in+1)).y + (*(in+5)).z); 
     (*(out+1)).z -= kappa*((*(in+1)).z - (*(in+3)).y);
     (*(out+1)).w -= kappa*((*(in+1)).w + (*(in+3)).x); 
     
     (*(out+2)).x -= kappa*((*(in+2)).x - (*(in+3)).w);
     (*(out+2)).y -= kappa*((*(in+2)).y + (*(in+3)).z);
     (*(out+2)).z -= kappa*((*(in+2)).z - (*(in+4)).y);
     (*(out+2)).w -= kappa*((*(in+2)).w + (*(in+4)).x);     
     
     (*(out+3)).x -= kappa*((*(in+3)).x + (*(in+1)).w);
     (*(out+3)).y -= kappa*((*(in+3)).y - (*(in+1)).z);     
     (*(out+3)).z -= kappa*((*(in+3)).z + (*(in+2)).y);
     (*(out+3)).w -= kappa*((*(in+3)).w - (*(in+2)).x);       
     
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+0)).y);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+0)).x);    
     (*(out+4)).x -= kappa*((*(in+4)).x + (*(in+2)).w);
     (*(out+4)).y -= kappa*((*(in+4)).y - (*(in+2)).z);     

     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+0)).w);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+0)).z);     
     (*(out+5)).z -= kappa*((*(in+5)).z + (*(in+1)).y);
     (*(out+5)).w -= kappa*((*(in+5)).w - (*(in+1)).x);     
     
}


//-kappa(r + gamma_mu)
__device__ void dev_kappaP1_minus(dev_spinor * out, dev_spinor * in, REAL kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+4)).w);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+4)).z);
     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+5)).y);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+5)).x);    

     (*(out+1)).x -= kappa*((*(in+1)).x + (*(in+5)).w);
     (*(out+1)).y -= kappa*((*(in+1)).y - (*(in+5)).z); 
     (*(out+1)).z -= kappa*((*(in+1)).z + (*(in+3)).y);
     (*(out+1)).w -= kappa*((*(in+1)).w - (*(in+3)).x); 
     
     (*(out+2)).x -= kappa*((*(in+2)).x + (*(in+3)).w);
     (*(out+2)).y -= kappa*((*(in+2)).y - (*(in+3)).z);
     (*(out+2)).z -= kappa*((*(in+2)).z + (*(in+4)).y);
     (*(out+2)).w -= kappa*((*(in+2)).w - (*(in+4)).x);     
     
     (*(out+3)).x -= kappa*((*(in+3)).x - (*(in+1)).w);
     (*(out+3)).y -= kappa*((*(in+3)).y + (*(in+1)).z);     
     (*(out+3)).z -= kappa*((*(in+3)).z - (*(in+2)).y);
     (*(out+3)).w -= kappa*((*(in+3)).w + (*(in+2)).x);       
     
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+0)).y);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+0)).x);    
     (*(out+4)).x -= kappa*((*(in+4)).x - (*(in+2)).w);
     (*(out+4)).y -= kappa*((*(in+4)).y + (*(in+2)).z);     

     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+0)).w);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+0)).z);     
     (*(out+5)).z -= kappa*((*(in+5)).z - (*(in+1)).y);
     (*(out+5)).w -= kappa*((*(in+5)).w + (*(in+1)).x);     
     
}





//-kappa(r - gamma_mu)
__device__ void dev_kappaP2_plus(dev_spinor * out, dev_spinor * in, REAL kappa){


     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+4)).z);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+4)).w);
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+0)).x);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+0)).y);    
     
 
     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+5)).x);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+5)).y);
     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+0)).z);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+0)).w);     
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x + (*(in+5)).z);
     (*(out+1)).y -= kappa*( (*(in+1)).y + (*(in+5)).w);
     (*(out+5)).z -= kappa*( (*(in+5)).z + (*(in+1)).x);
     (*(out+5)).w -= kappa*( (*(in+5)).w + (*(in+1)).y);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z - (*(in+3)).x);
     (*(out+1)).w -= kappa*( (*(in+1)).w - (*(in+3)).y);
     (*(out+3)).x -= kappa*( (*(in+3)).x - (*(in+1)).z);
     (*(out+3)).y -= kappa*( (*(in+3)).y - (*(in+1)).w);     
     
     
     (*(out+2)).x -= kappa*( (*(in+2)).x - (*(in+3)).z);
     (*(out+2)).y -= kappa*( (*(in+2)).y - (*(in+3)).w);
     (*(out+3)).z -= kappa*( (*(in+3)).z - (*(in+2)).x);
     (*(out+3)).w -= kappa*( (*(in+3)).w - (*(in+2)).y);     
     
     
     (*(out+2)).z -= kappa*( (*(in+2)).z - (*(in+4)).x);
     (*(out+2)).w -= kappa*( (*(in+2)).w - (*(in+4)).y);
     (*(out+4)).x -= kappa*( (*(in+4)).x - (*(in+2)).z);
     (*(out+4)).y -= kappa*( (*(in+4)).y - (*(in+2)).w);
   
     
}


//-kappa(r + gamma_mu)  kappa reell !!!!
__device__ void dev_kappaP2_minus(dev_spinor * out, dev_spinor * in, REAL kappa){


     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+4)).z);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+4)).w);
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+0)).x);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+0)).y);    
     
 
     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+5)).x);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+5)).y);
     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+0)).z);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+0)).w);     
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x - (*(in+5)).z);
     (*(out+1)).y -= kappa*( (*(in+1)).y - (*(in+5)).w);
     (*(out+5)).z -= kappa*( (*(in+5)).z - (*(in+1)).x);
     (*(out+5)).w -= kappa*( (*(in+5)).w - (*(in+1)).y);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z + (*(in+3)).x);
     (*(out+1)).w -= kappa*( (*(in+1)).w + (*(in+3)).y);
     (*(out+3)).x -= kappa*( (*(in+3)).x + (*(in+1)).z);
     (*(out+3)).y -= kappa*( (*(in+3)).y + (*(in+1)).w);     
     
     
     (*(out+2)).x -= kappa*( (*(in+2)).x + (*(in+3)).z);
     (*(out+2)).y -= kappa*( (*(in+2)).y + (*(in+3)).w);
     (*(out+3)).z -= kappa*( (*(in+3)).z + (*(in+2)).x);
     (*(out+3)).w -= kappa*( (*(in+3)).w + (*(in+2)).y);     
     
     
     (*(out+2)).z -= kappa*( (*(in+2)).z + (*(in+4)).x);
     (*(out+2)).w -= kappa*( (*(in+2)).w + (*(in+4)).y);
     (*(out+4)).x -= kappa*( (*(in+4)).x + (*(in+2)).z);
     (*(out+4)).y -= kappa*( (*(in+4)).y + (*(in+2)).w);
   
     
}



//-kappa(r - gamma_mu) kappa reell !!!!
__device__ void dev_kappaP3_plus(dev_spinor * out, dev_spinor * in, REAL kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x - (*(in+3)).y);
     (*(out+0)).y -= kappa*( (*(in+0)).y + (*(in+3)).x);
     (*(out+3)).x -= kappa*( (*(in+3)).x + (*(in+0)).y);
     (*(out+3)).y -= kappa*( (*(in+3)).y - (*(in+0)).x);    
     

     (*(out+0)).z -= kappa*( (*(in+0)).z - (*(in+3)).w);
     (*(out+0)).w -= kappa*( (*(in+0)).w + (*(in+3)).z);
     (*(out+3)).z -= kappa*( (*(in+3)).z + (*(in+0)).w);
     (*(out+3)).w -= kappa*( (*(in+3)).w - (*(in+0)).z);    
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x - (*(in+4)).y);
     (*(out+1)).y -= kappa*( (*(in+1)).y + (*(in+4)).x);
     (*(out+4)).x -= kappa*( (*(in+4)).x + (*(in+1)).y);
     (*(out+4)).y -= kappa*( (*(in+4)).y - (*(in+1)).x);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z + (*(in+4)).w);
     (*(out+1)).w -= kappa*( (*(in+1)).w - (*(in+4)).z);
     (*(out+4)).z -= kappa*( (*(in+4)).z - (*(in+1)).w);
     (*(out+4)).w -= kappa*( (*(in+4)).w + (*(in+1)).z);     
     
       
     (*(out+2)).x -= kappa*( (*(in+2)).x + (*(in+5)).y);
     (*(out+2)).y -= kappa*( (*(in+2)).y - (*(in+5)).x);
     (*(out+5)).x -= kappa*( (*(in+5)).x - (*(in+2)).y);
     (*(out+5)).y -= kappa*( (*(in+5)).y + (*(in+2)).x);    
     

     (*(out+2)).z -= kappa*( (*(in+2)).z + (*(in+5)).w);
     (*(out+2)).w -= kappa*( (*(in+2)).w - (*(in+5)).z);
     (*(out+5)).z -= kappa*( (*(in+5)).z - (*(in+2)).w);
     (*(out+5)).w -= kappa*( (*(in+5)).w + (*(in+2)).z);
  
}


//-kappa(r + gamma_mu) kappa reell !!!
__device__ void dev_kappaP3_minus(dev_spinor * out, dev_spinor * in, REAL kappa){

     (*(out+0)).x -= kappa*( (*(in+0)).x + (*(in+3)).y);
     (*(out+0)).y -= kappa*( (*(in+0)).y - (*(in+3)).x);
     (*(out+3)).x -= kappa*( (*(in+3)).x - (*(in+0)).y);
     (*(out+3)).y -= kappa*( (*(in+3)).y + (*(in+0)).x);    
     

     (*(out+0)).z -= kappa*( (*(in+0)).z + (*(in+3)).w);
     (*(out+0)).w -= kappa*( (*(in+0)).w - (*(in+3)).z);
     (*(out+3)).z -= kappa*( (*(in+3)).z - (*(in+0)).w);
     (*(out+3)).w -= kappa*( (*(in+3)).w + (*(in+0)).z);    
     
     
     (*(out+1)).x -= kappa*( (*(in+1)).x + (*(in+4)).y);
     (*(out+1)).y -= kappa*( (*(in+1)).y - (*(in+4)).x);
     (*(out+4)).x -= kappa*( (*(in+4)).x - (*(in+1)).y);
     (*(out+4)).y -= kappa*( (*(in+4)).y + (*(in+1)).x);     
     
     
     (*(out+1)).z -= kappa*( (*(in+1)).z - (*(in+4)).w);
     (*(out+1)).w -= kappa*( (*(in+1)).w + (*(in+4)).z);
     (*(out+4)).z -= kappa*( (*(in+4)).z + (*(in+1)).w);
     (*(out+4)).w -= kappa*( (*(in+4)).w - (*(in+1)).z);     
     
       
     (*(out+2)).x -= kappa*( (*(in+2)).x - (*(in+5)).y);
     (*(out+2)).y -= kappa*( (*(in+2)).y + (*(in+5)).x);
     (*(out+5)).x -= kappa*( (*(in+5)).x + (*(in+2)).y);
     (*(out+5)).y -= kappa*( (*(in+5)).y - (*(in+2)).x);    
     

     (*(out+2)).z -= kappa*( (*(in+2)).z - (*(in+5)).w);
     (*(out+2)).w -= kappa*( (*(in+2)).w + (*(in+5)).z);
     (*(out+5)).z -= kappa*( (*(in+5)).z + (*(in+2)).w);
     (*(out+5)).w -= kappa*( (*(in+5)).w - (*(in+2)).z);
  
}







//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_plus(dev_spinor * out, dev_spinor * in, dev_complex kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
  
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
        
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}






//-kappa(r - gamma_mu)
__device__ void dev_kappaP0_minus(dev_spinor * out, dev_spinor * in, dev_complex kappa){


     (*(out+0)).x -= (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+0)).y -= (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im;
     (*(out+0)).x += (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+0)).y += (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;
     
     (*(out+3)).x -= (*(in+3)).x*kappa.re - (*(in+3)).y*kappa.im;
     (*(out+3)).y -= (*(in+3)).y*kappa.re + (*(in+3)).x*kappa.im;   
     (*(out+3)).x += (*(in+0)).x*kappa.re - (*(in+0)).y*kappa.im;
     (*(out+3)).y += (*(in+0)).y*kappa.re + (*(in+0)).x*kappa.im; 



     (*(out+0)).z -= (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+0)).w -= (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;
     (*(out+0)).z += (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+0)).w += (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;     
          
     (*(out+3)).z -= (*(in+3)).z*kappa.re - (*(in+3)).w*kappa.im;
     (*(out+3)).w -= (*(in+3)).w*kappa.re + (*(in+3)).z*kappa.im;
     (*(out+3)).z += (*(in+0)).z*kappa.re - (*(in+0)).w*kappa.im;
     (*(out+3)).w += (*(in+0)).w*kappa.re + (*(in+0)).z*kappa.im;

 
 
     (*(out+1)).x -= (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+1)).y -= (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im;
     (*(out+1)).x += (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+1)).y += (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;     
     
     (*(out+4)).x -= (*(in+4)).x*kappa.re - (*(in+4)).y*kappa.im;
     (*(out+4)).y -= (*(in+4)).y*kappa.re + (*(in+4)).x*kappa.im;    
     (*(out+4)).x += (*(in+1)).x*kappa.re - (*(in+1)).y*kappa.im;
     (*(out+4)).y += (*(in+1)).y*kappa.re + (*(in+1)).x*kappa.im; 
 
 
 
     (*(out+1)).z -= (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+1)).w -= (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;
     (*(out+1)).z += (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+1)).w += (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;     
     
     (*(out+4)).z -= (*(in+4)).z*kappa.re - (*(in+4)).w*kappa.im;
     (*(out+4)).w -= (*(in+4)).w*kappa.re + (*(in+4)).z*kappa.im;    
     (*(out+4)).z += (*(in+1)).z*kappa.re - (*(in+1)).w*kappa.im;
     (*(out+4)).w += (*(in+1)).w*kappa.re + (*(in+1)).z*kappa.im;      
 

     
     (*(out+2)).x -= (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+2)).y -= (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im;
     (*(out+2)).x += (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+2)).y += (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;
                
     (*(out+5)).x -= (*(in+5)).x*kappa.re - (*(in+5)).y*kappa.im;
     (*(out+5)).y -= (*(in+5)).y*kappa.re + (*(in+5)).x*kappa.im;   
     (*(out+5)).x += (*(in+2)).x*kappa.re - (*(in+2)).y*kappa.im;
     (*(out+5)).y += (*(in+2)).y*kappa.re + (*(in+2)).x*kappa.im; 
   
   
       
     (*(out+2)).z -= (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+2)).w -= (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im;
     (*(out+2)).z += (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+2)).w += (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;
               
     (*(out+5)).z -= (*(in+5)).z*kappa.re - (*(in+5)).w*kappa.im;
     (*(out+5)).w -= (*(in+5)).w*kappa.re + (*(in+5)).z*kappa.im;   
     (*(out+5)).z += (*(in+2)).z*kappa.re - (*(in+2)).w*kappa.im;
     (*(out+5)).w += (*(in+2)).w*kappa.re + (*(in+2)).z*kappa.im; 
  
}











//applies the Hopping Part Even-Odd !
//the gauge field is the complete gaugefield!
//the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
//from 0..VOLUME/2-1, eo = 0 or 1
//the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
//the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
//other eo-sublattice.
//for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
//the boundary conditions are implemented as in Hopping_Matrix.c
//mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)  
__global__ void dev_Hopping_Matrix(const dev_su3_2v * gf, const dev_spinor * sin, dev_spinor * sout, const int * gfindex_site, const int* gfindex_nextsite, const int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3_pad gfsmem[BLOCK];



  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  
  
  if(pos < dev_VOLUME){
  

  dev_zero_spinor(&(ssum[0])); // zero sum        
 #ifdef TEMPORALGAUGE
  int spatialvol = dev_LX*dev_LY*dev_LZ;
 #endif


//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
             //hoppos = tex1Dfetch(nn_tex,8*pos);
            //color
            
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              
              if((gfindex_site[pos]/spatialvol) != (dev_T-1) ){
              #ifdef USETEXTURE
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
              #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
              #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos]),&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif
              }
            #else
              #ifdef GF_8
              dev_reconstructgf_8texref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref(gf, 4*(gfindex_site[pos]),&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
              #endif
            #endif
            
            //-kappa(r - gamma_mu)
            #ifdef GF_8
             dev_kappaP0_plus(&(ssum[0]), &(shelp1[0]), dev_cconj(dev_k0));
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
              dev_Gamma0(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));
	    #endif
	    
//l==0,t
            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+4);
            //color
            #ifdef TEMPORALGAUGE
              // gf == ID for t != T-1 => just read the spinor
              if((gfindex_nextsite[hoppos]/spatialvol) != (dev_T-1) ){
               #ifdef USETEXTURE
                shelp1[0] = tex1Dfetch(spin_tex,6*hoppos);
                shelp1[1] = tex1Dfetch(spin_tex,6*hoppos+1);
                shelp1[2] = tex1Dfetch(spin_tex,6*hoppos+2);
                shelp1[3] = tex1Dfetch(spin_tex,6*hoppos+3);
                shelp1[4] = tex1Dfetch(spin_tex,6*hoppos+4);
                shelp1[5] = tex1Dfetch(spin_tex,6*hoppos+5);
               #else
                shelp1[0] = sin[6*hoppos];
                shelp1[1] = sin[6*hoppos+1];
                shelp1[2] = sin[6*hoppos+2];
                shelp1[3] = sin[6*hoppos+3];
                shelp1[4] = sin[6*hoppos+4];
                shelp1[5] = sin[6*hoppos+5];
               #endif
              }
              else{
                // gf != ID for t == T-1 => mult spinor with gf
                #ifdef GF_8
                dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #else
                dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
                #endif
                #ifdef USETEXTURE
                  dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
                #else
                  dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
                #endif 
              }
            #else            
              #ifdef GF_8
              dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #else
              dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos],&(gfsmem[ix].m));
              #endif
              #ifdef USETEXTURE
                dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));  
              #else
                dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
              #endif 
            #endif
            
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP0_minus(&(ssum[0]), &(shelp1[0]), dev_k0);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
              dev_Gamma0(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            #endif




//l==3,z 
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
             //hoppos = tex1Dfetch(nn_tex,8*pos+3);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf, 4*(gfindex_site[pos])+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)    
            #ifdef GF_8
            dev_kappaP3_plus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));
	    #endif
//l==3,z               
            
            //negative direction
            hoppos = nn_evenodd[8*pos+7];
             //hoppos = tex1Dfetch(nn_tex,8*pos+7); 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP3_minus(&(ssum[0]), &(shelp1[0]), dev_k3.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
              dev_Gamma3(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            #endif




//l==2,y 
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
             //hoppos = tex1Dfetch(nn_tex,8*pos+2);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos])+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP2_plus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            #endif

//l==2,y        

            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+6);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP2_minus(&(ssum[0]), &(shelp1[0]), dev_k2.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
              dev_Gamma2(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
	    #endif



//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
             //hoppos = tex1Dfetch(nn_tex,8*pos+1);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref(gf,4*(gfindex_site[pos])+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            #ifdef GF_8
              dev_kappaP1_plus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));
	    #endif


//l==1,x 
            
            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
             //hoppos = tex1Dfetch(nn_tex,8*pos+5);
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix].m));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix].m, hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix].m, &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            #ifdef GF_8
              dev_kappaP1_minus(&(ssum[0]), &(shelp1[0]), dev_k1.re);
            #else
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
              dev_Gamma1(&(shelp1[0]));
              dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));      
            #endif
 
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos])); 
  }
}



/*

//applies the Hopping Part Even-Odd !
//the gauge field is the complete gaugefield!
//the gauge field at the local point is reconstructed by 2*pos+eo where pos is the eo-position
//from 0..VOLUME/2-1, eo = 0 or 1
//the positions in the gauge fields are passed in "gfindex_site" for gf's that are attached at
//the actual positions and in "gfindex_nextsite" for gf's that start at a position of the 
//other eo-sublattice.
//for the hopping positions of the eo-spinor field we use on of the two dedicated eo-nn fields
//the boundary conditions are implemented as in Hopping_Matrix.c
//mult with complex conjugate k0,k1,k2,k3 in positive direction because
// psi(x+mu) != exp(i theta_mu) psi(x)  
__global__ void dev_Hopping_Matrix(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * gfindex_site,int* gfindex_nextsite, int * nn_evenodd, const int eo){

  int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){

  dev_zero_spinor(&(ssum[0])); // zero sum        
//hopping term                
//l==0,t
            //positive direction
            hoppos = nn_evenodd[8*pos];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos]),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos]),&(gfsmem[ix]));
            #endif
            
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+4]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos],&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos],&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));     
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = nn_evenodd[8*pos+3];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(3),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+7]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(3),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = nn_evenodd[8*pos+2];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(2),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = nn_evenodd[8*pos+6]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(2),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = nn_evenodd[8*pos+1];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(4*(gfindex_site[pos])+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(4*(gfindex_site[pos])+(1),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r - gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = nn_evenodd[8*pos+5]; 
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(4*gfindex_nextsite[hoppos]+(1),&(gfsmem[ix]));
            #endif
            dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            //-kappa(r + gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));      
               
 
        //copy to output spinor
        dev_copy_spinor(&(ssum[0]),&(sout[6*pos])); 
  }
}



*/



// aequivalent to Qtm_pm_psi in tm_operators.c
extern "C" void dev_Qtm_pm_psi(dev_spinor* spinin, dev_spinor* spinout, int gridsize, int blocksize, int gridsize2, int blocksize2){
  //spinin == odd
  //spinout == odd
  
  //Q_{-}
  #ifdef USETEXTURE
    bind_texture_spin(spinin,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0           
  //unbind_texture_nn();           
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,dev_spin_eo2, -1.);
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
            (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(spinin, dev_spin_eo1,  dev_spin_eo2, -1.);
  
  //Q_{+}
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2,1);
  #endif
  //bind_texture_nn(dev_nn_eo);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
          (dev_gf, dev_spin_eo2, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
  //unbind_texture_nn();      
  #ifdef USETEXTURE  
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_inv<<<gridsize2, blocksize2>>>(dev_spin_eo1,spinout, +1.);
  
  #ifdef USETEXTURE
    bind_texture_spin(spinout,1);
  #endif
  //bind_texture_nn(dev_nn_oe);
  //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<gridsize, blocksize>>>
             (dev_gf, spinout, dev_spin_eo1, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();  
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  dev_mul_one_pm_imu_sub_mul_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo2, dev_spin_eo1,  spinout , +1.); 
}






//applies the full tm Operator
// uses texture cache (spin_tex) for input spinor
// runs through whole lattice for output spinor
// D_psi uses phase_mu and not ka_mu for the boundary conds (vice versa in HoppingMatrix) 
// -> thats why complexmult and complexcgmult are interchanged in dev_HoppingMatrix and in 
// dev_tm_dirac_kappa
__global__ void dev_tm_dirac_kappa(dev_su3_2v * gf, dev_spinor * sin, dev_spinor * sout, int * dev_nn){
    int pos,hoppos;
    dev_spinor shelp1[6], ssum[6];
    __shared__ dev_su3 gfsmem[BLOCK];
    

  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  int ix = threadIdx.x;
  if(pos < dev_VOLUME){
        
          //dev_zero_spinor(&(ssum[0])); // zero sum
          //skalarer Term
         #ifdef USETEXTURE
          ssum[0] = tex1Dfetch(spin_tex,6*pos);
          ssum[1] = tex1Dfetch(spin_tex,6*pos+1);
          ssum[2] = tex1Dfetch(spin_tex,6*pos+2);
          ssum[3] = tex1Dfetch(spin_tex,6*pos+3);
          ssum[4] = tex1Dfetch(spin_tex,6*pos+4);
          ssum[5] = tex1Dfetch(spin_tex,6*pos+5);
	 #else
	  ssum[0] = sin[6*pos];
          ssum[1] = sin[6*pos+1];
          ssum[2] = sin[6*pos+2];
          ssum[3] = sin[6*pos+3];
          ssum[4] = sin[6*pos+4];
          ssum[5] = sin[6*pos+5];
	 #endif
          
//hopping term                
//l==0,t
            //positive direction
            hoppos = dev_nn[8*pos];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*pos,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,4*pos,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k0,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+4];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*hoppos,&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*hoppos,&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));  
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif    
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(0,&(shelp1[0]));
            dev_Gamma0(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk0,&(shelp1[0]), &(ssum[0]));


//l==3,z               
            //positive direction
            hoppos = dev_nn[8*pos+3];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*pos+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,4*pos+(3),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k3,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+7];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*hoppos+(3),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*hoppos+(3),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(3,&(shelp1[0]));
            dev_Gamma3(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk3,&(shelp1[0]), &(ssum[0]));
         
         
//l==2,y        
            //positive direction
            hoppos = dev_nn[8*pos+2];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*pos+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,4*pos+(2),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k2,&(shelp1[0]), &(ssum[0]));
            
            //negative direction
            hoppos = dev_nn[8*pos+6];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*hoppos+(2),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*hoppos+(2),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(2,&(shelp1[0]));
            dev_Gamma2(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk2,&(shelp1[0]), &(ssum[0]));


//l==1,x 
            //positive direction
            hoppos = dev_nn[8*pos+1];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref(gf,4*pos+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref(gf,4*pos+(1),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r - gamma_mu)
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexmult_add_assign_spinor(&(ssum[0]),dev_k1,&(shelp1[0]), &(ssum[0]));

            //negative direction
            hoppos = dev_nn[8*pos+5];
            //color
            #ifdef GF_8
            dev_reconstructgf_8texref_dagger(gf,4*hoppos+(1),&(gfsmem[ix]));
            #else
            dev_reconstructgf_2vtexref_dagger(gf,4*hoppos+(1),&(gfsmem[ix]));
            #endif
            #ifdef USETEXTURE
              dev_su3MtV_spintex(gfsmem[ix], hoppos, &(shelp1[0]));
            #else
              dev_su3MtV(gfsmem[ix], &(sin[6*hoppos]), &(shelp1[0]));
            #endif
            //-kappa(r + gamma_mu)
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));
            //dev_GammatV(1,&(shelp1[0]));
            dev_Gamma1(&(shelp1[0]));
            dev_complexcgmult_add_assign_spinor(&(ssum[0]),dev_mk1,&(shelp1[0]), &(ssum[0]));  
          
          
          
          //gamma5 term
         #ifdef USETEXTURE
          shelp1[0] = tex1Dfetch(spin_tex,6*pos);
          shelp1[1] = tex1Dfetch(spin_tex,6*pos+1);
          shelp1[2] = tex1Dfetch(spin_tex,6*pos+2);
          shelp1[3] = tex1Dfetch(spin_tex,6*pos+3);
          shelp1[4] = tex1Dfetch(spin_tex,6*pos+4);
          shelp1[5] = tex1Dfetch(spin_tex,6*pos+5);
         #else
          shelp1[0] = sin[6*pos];
          shelp1[1] = sin[6*pos+1];
          shelp1[2] = sin[6*pos+2];
          shelp1[3] = sin[6*pos+3];
          shelp1[4] = sin[6*pos+4];
          shelp1[5] = sin[6*pos+5];
         #endif 
          
          
          //dev_GammatV(4,&(shelp1[0]));
          dev_Gamma5(&(shelp1[0]));
          dev_complexmult_add_assign_spinor(&(ssum[0]),dev_initcomplex(0.0,2.0*kappa*mu),&(shelp1[0]), &(sout[6*pos]));
  }
}






extern "C" void dev_tm_dirac_dagger_kappa(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
 int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){
 int gridsize;
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 dim3 blockdim2(128,1,1);
 dim3 blockdim(xsize,ysize);
 
  dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int)(VOLUME/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
  dev_gamma5 <<<griddim2, blockdim2 >>> (spinin,spinout);
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spinout, spinin, dev_nn);
  dev_gamma5 <<<griddim2, blockdim2 >>>(spinin,spinout);
}







__device__ inline REAL dev_skalarprod_spinor(dev_spinor * s1, dev_spinor * s2){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s2+i)).x + (*(s1+i)).y*(*(s2+i)).y + (*(s1+i)).z*(*(s2+i)).z + (*(s1+i)).w*(*(s2+i)).w);
  }
  return skalprod;
}




__device__ inline REAL dev_squarenorm_spinor(dev_spinor * s1){
  REAL skalprod = 0.0;
  int i;
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    skalprod += ((*(s1+i)).x*(*(s1+i)).x + (*(s1+i)).y*(*(s1+i)).y + (*(s1+i)).z*(*(s1+i)).z + (*(s1+i)).w*(*(s1+i)).w);
  }
  return skalprod;
}



__device__ inline REAL dev_squarenorm_spinor_tex(int pos){
  REAL skalprod = 0.0;
  int i;
  float4 help;
  
  #pragma unroll 6
  for(i=0;i<6;i++){ //color + spin
    help = tex1Dfetch(spin_tex2,6*pos+i);
    skalprod += help.x*help.x + help.y*help.y + help.z*help.z + help.w*help.w;
  }
  return skalprod;
}




//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_skalarprod_spinor_field2(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}








//only 1 dim parallel possible, because need __syncthread !
__global__ void dev_squarenorm_spinor_field(dev_spinor* s1, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth;
  REAL ks,kc,ds,tr,ts,tt;
  
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(blockDim.x > dev_VOLUME){
     stepwidth = 1;  
   }
   else{
     stepwidth = dev_VOLUME/(gridDim.x*blockDim.x);
   }
   
     int start = (blockIdx.x*blockDim.x + threadIdx.x)*stepwidth;
     int end = (blockIdx.x*blockDim.x + threadIdx.x+1)*stepwidth;  
   
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
         //ds = dev_squarenorm_spinor_tex(pos);
         ds = dev_squarenorm_spinor(&(s1[6*pos]));
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x] = kc;
    __syncthreads();
    
    
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (erg[blockIdx.x])=kc;
    }//threadIdx==0;
}






//only 1 dim parallel, because need __syncthread !
__global__ void dev_skalarprod_spinor_field(dev_spinor* s1, dev_spinor* s2, REAL* erg){
  __shared__ REAL shrinkarray[ACCUM_N];
  int pos,stepwidth, sweepsperthread;
  REAL ks,kc,ds,tr,ts,tt;
   
   // ADD ERROR HERE if t > maxblockdim
   
   ks=0.0;
   kc=0.0; 
   
   if(ACCUM_N > dev_VOLUME){
     stepwidth = 1;
     sweepsperthread = 1;  
   }
   else{
     stepwidth = dev_VOLUME/ACCUM_N;
     sweepsperthread = ACCUM_N/blockDim.x;
   }
    
   
   
 for(int j = 0; j < sweepsperthread; j++){
   
     int start = (threadIdx.x + j*blockDim.x)*stepwidth;
     int end = (threadIdx.x+j*blockDim.x+1)*stepwidth;  
     ks=0.0;
     kc=0.0; 
     
   for(pos=start;pos<end; pos++){
     if(pos < dev_VOLUME){
          ds = dev_skalarprod_spinor(&(s1[6*pos]), &(s2[6*pos] ));
          
      }
      else{
        ds = 0.0;
      } 
          // Kahan summation
          tr=ds+kc;
          ts=tr+ks;
          tt=ts-ks;
          ks=ts;
          kc=tr-tt;
    }
    kc=ks+kc;
    shrinkarray[threadIdx.x+j*blockDim.x] = kc;
  }
  __syncthreads();
   
 
    for(int stride = ACCUM_N / 2; stride > 0; stride >>= 1){
       __syncthreads();
       for(int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x)
           shrinkarray[iAccum] += shrinkarray[stride + iAccum];
    }

    if(threadIdx.x == 0) (*erg) = shrinkarray[0];
    
    
    /*
    if(threadIdx.x==0){
      ks=0.0;
      kc=0.0; 
      int k;
      for(k=0; k<sweepsperthread*blockDim.x; k++){
        ds = shrinkarray[k];
        tr=ds+kc;
        ts=tr+ks;
        tt=ts-ks;
        ks=ts;
        kc=tr-tt;
      }
      kc=ks+kc;
      (*erg)=kc;
    }//threadIdx==0;
    */
    
}




__global__ void dev_zero_spinor_field(dev_spinor* s1){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
          dev_zero_spinor(&(s1[6*pos]));
  }
}




__global__ void dev_copy_spinor_field(dev_spinor* s1, dev_spinor* s2){
    int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
      dev_copy_spinor(&(s1[6*pos]),&(s2[6*pos]));
  } 
}



__global__ void dev_skalarmult_add_assign_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* s2, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_add_assign_spinor(&(s1[6*pos]), lambda ,&(s2[6*pos]), &(so[6*pos]) );
  }
}



__global__ void dev_skalarmult_spinor_field(dev_spinor* s1, REAL lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[6*pos]), dev_initcomplex(lambda,0.0) , &(so[6*pos]) );
  }
}  



__global__ void dev_complexmult_spinor_field(dev_spinor* s1, dev_complex lambda, dev_spinor* so){
  int pos;
  pos= threadIdx.x + blockDim.x*blockIdx.x;  
  if(pos < dev_VOLUME){
    dev_skalarmult_spinor(&(s1[6*pos]), lambda , &(so[6*pos]) );
  }
}






// init the gpu inner solver, assigen constants etc.
__global__ void he_cg_init (int* grid, REAL param_kappa, REAL param_mu, dev_complex k0, dev_complex k1, dev_complex k2, dev_complex k3){
  dev_LX = grid[0];
  dev_LY = grid[1];
  dev_LZ = grid[2];
  dev_T = grid[3];
  dev_VOLUME = grid[4]; // grid[4] is initialized 1/2 VOLUME for eo
  
  kappa = param_kappa;
  mu = param_mu;
  twokappamu = 2.0*param_kappa*param_mu;
  
  dev_k0.re = k0.re;
  dev_k0.im = k0.im;
  dev_mk0.re = -k0.re;
  dev_mk0.im = -k0.im;
  
  dev_k1.re = k1.re;
  dev_k1.im = k1.im;
  dev_mk1.re = -k1.re;
  dev_mk1.im = -k1.im;
  
  dev_k2.re = k2.re;
  dev_k2.im = k2.im;
  dev_mk2.re = -k2.re;
  dev_mk2.im = -k2.im;
  
  dev_k3.re = k3.re;
  dev_k3.im = k3.im;
  dev_mk3.re = -k3.re;
  dev_mk3.im = -k3.im;
}




// code to list available devices, not yet included in main code
// this is copied from the CUDA sdk 
extern "C" int find_devices(){
int deviceCount, dev;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
        printf("There is no device supporting CUDA\n");
    for (dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
        }
        printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
        printf("  Major revision number:                         %d\n",
               deviceProp.major);
        printf("  Minor revision number:                         %d\n",
               deviceProp.minor);
        printf("  Total amount of global memory:                 %u bytes\n",
               deviceProp.totalGlobalMem);
    #if CUDART_VERSION >= 2000
        printf("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
        printf("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
    #endif
        printf("  Total amount of constant memory:               %u bytes\n",
               deviceProp.totalConstMem); 
        printf("  Total amount of shared memory per block:       %u bytes\n",
               deviceProp.sharedMemPerBlock);
        printf("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
        printf("  Warp size:                                     %d\n",
               deviceProp.warpSize);
        printf("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf("  Maximum memory pitch:                          %u bytes\n",
               deviceProp.memPitch);
        printf("  Texture alignment:                             %u bytes\n",
               deviceProp.textureAlignment);
        printf("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
    #if CUDART_VERSION >= 2000
        printf("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
    #endif
    }
    return(deviceCount);
}






extern "C" int bind_texture_spin(dev_spinor* s, int i){
  
  size_t size;
  if(even_odd_flag){
    size = sizeof(float4)*6*VOLUME/2;
  }
  else{
    size = sizeof(float4)*6*VOLUME;
  }
   
  
  switch(i){
    case 1:
      //printf("Binding texture to spinorfield 1\n");
      spin_texRefPtr = NULL;
      cudaGetTextureReference(&spin_texRefPtr, "spin_tex");
      spin_channelDesc =  cudaCreateChannelDesc<float4>();
      cudaBindTexture(0, spin_texRefPtr, s, &spin_channelDesc, size);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
      return(0);
   
    case 2:
      //printf("Binding texture to spinorfield 2\n");
      spin_texRefPtr2 = NULL;
      cudaGetTextureReference(&spin_texRefPtr2, "spin_tex2");
      spin_channelDesc2 =  cudaCreateChannelDesc<float4>();
      cudaBindTexture(0, spin_texRefPtr2, s, &spin_channelDesc2, size);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));  
      return(0);
  }
return(1);  
}


extern "C" int unbind_texture_spin(int i){
  switch(i){
    case 1:
      //printf("Unbinding texture of spinorfield 1\n");
      cudaUnbindTexture(spin_texRefPtr);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
      return(0);
    case 2:
      //printf("Unbinding texture of spinorfield 2\n");
      cudaUnbindTexture(spin_texRefPtr2);
      //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
      return(0);    
  }
  
return(1);
}



extern "C" int bind_texture_gf(dev_su3_2v * gf){
 //printf("Binding texture to gaugefield\n");
 
 #ifdef GF_8
 size_t size = sizeof(float4)*2*VOLUME*4;
 #else
 size_t size = sizeof(float4)*3*VOLUME*4;
 #endif
 
 cudaGetTextureReference(&gf_texRefPtr, "gf_tex");
 gf_channelDesc =  cudaCreateChannelDesc<float4>();
 cudaBindTexture(0, gf_texRefPtr, gf, &gf_channelDesc, size);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_gf(){
 //printf("Unbinding texture to gaugefield\n");
 cudaUnbindTexture(gf_texRefPtr);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}



extern "C" int bind_texture_nn(int* nn){
 //printf("Binding texture to nn field\n");
  size_t size;
  if(even_odd_flag){
    size = sizeof(int)*8*VOLUME/2;
  }
  else{
    size = sizeof(int)*8*VOLUME;
  }
 

 cudaGetTextureReference(&nn_texRefPtr, "nn_tex");
 nn_channelDesc =  cudaCreateChannelDesc<int>();
 cudaBindTexture(0, nn_texRefPtr, nn, &nn_channelDesc, size);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}


extern "C" int unbind_texture_nn(){
 //printf("Unbinding texture to nn field\n");
 cudaUnbindTexture(nn_texRefPtr);
 //printf("%s\n", cudaGetErrorString(cudaGetLastError()));    
 return(0);
}






extern "C" void test_operator(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize){
 
 int  gridsize;

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 
 dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int) VOLUME/BLOCK + 1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
 
 
  dev_complex h0,h1,h2,h3;
  h0.re = (REAL)ka0.re;    h0.im = (REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = (REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = (REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = (REAL)ka3.im;
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
 
 
  REAL scaleparam = sqrt(1.0/(2.0 * (REAL) hostkappa));
  dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam*scaleparam, spin4);
 
 #ifdef USETEXTURE
   bind_texture_gf(gf);
   bind_texture_spin(spin4,1);
 #endif 
  // apply D_tm
  dev_tm_dirac_kappa <<<griddim3, blockdim3 >>>(gf, spin4, spinout, nn_grid);

 #ifdef USETEXTURE
  unbind_texture_gf();
  unbind_texture_spin(1);
 #endif
}






extern "C" void dev_cg(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize, int rescalekappa){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
 int N_recalcres = 10; // after N_recalcres iterations calculate r = A x_k - b
 
 
 // initialize grid and block, make sure VOLUME is a multiple of blocksize 
 if(VOLUME%DOTPROD_DIM != 0){
   printf("Error: VOLUME is not a multiple of DOTPROD_DIM. Aborting...\n");
   exit(100); 
 }

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME >= 128){
   gridsize =VOLUME/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);
 
 
 //if(VOLUME%BLOCK != 0){
 //  printf("Error: VOLUME is not a multiple of BLOCK. Aborting...\n");
 //  exit(100);
 //}
 dim3 blockdim3(BLOCK,1,1);
 if( VOLUME >= BLOCK){
   gridsize = (int) (VOLUME/BLOCK) +1;
 }
 else{
   gridsize=1;
 }
 dim3 griddim3(gridsize,1,1); 
 
 size_t size2 = sizeof(float4)*6*VOLUME;
 
 #ifdef USETEXTURE
   //Bind texture gf
   bind_texture_gf(gf);
  //Bind texture spinor to spin4 (D_tm is always applied to spin4)
  bind_texture_spin(spin4,1);
 #endif
 
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3;
  h0.re = (REAL)ka0.re;    h0.im = (REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = (REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = (REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = (REAL)ka3.im;
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(REAL));
 cudaMalloc((void **) &dotprod2, sizeof(REAL));
 cudaMalloc((void **) &rk, sizeof(REAL));
 cudaMalloc((void **) &alpha, sizeof(REAL));
 cudaMalloc((void **) &beta, sizeof(REAL));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 printf("have initialized cublas\n");
 
 
 // go over to kappa (if wanted)
 REAL scaleparam = sqrt(1.0/(2.0 * (REAL)hostkappa));
 printf("1/2kappa = %.8f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 
 
 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot (24*VOLUME, (const float *)spinin, 1, (const float *)spinin, 1);
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
  printf("Entering cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // D Ddagger    --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
  // mu -> -mu for twisted term
  // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
     // GAMMA5, mu -> -mu
     dev_gamma5 <<<griddim2, blockdim2 >>> (spin2,spin4);
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm 
     dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
  #ifdef USETEXTURE
   unbind_texture_spin(1);
  #endif
     //GAMMA5 mu -> -mu
     dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin4);
     dev_swapmu <<<1,1>>> ();
  #ifdef USETEXTURE
   bind_texture_spin(spin4,1);
  #endif
     //D_tm
     dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
  
  //Here we have used the output spinor (spinout) to temporarly take the field and to 
  //copy it to the texture field (spin4)!!

  
 //alpha
  host_dotprod = cublasSdot (24*VOLUME, (const float *) spin2, 1,
            (const float *) spin3, 1);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy (24*VOLUME,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  

 //x(k+1);
 cublasSaxpy (24*VOLUME, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);

 printf("%s\n", cudaGetErrorString(cudaGetLastError()));

  //Abbruch?
  host_dotprod = cublasSdot (24*VOLUME, (const float *) spin0, 1,(const float *) spin0, 1);
  
 if ((host_dotprod <= eps*sourcesquarenorm)){//error-limit erreicht
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 cublasSscal (24*VOLUME, host_beta, (float *)spin2, 1);
 cublasSaxpy (24*VOLUME, 1.0, (const float *) spin0,  1, (float *) spin2, 1);

 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
      //GAMMA5
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
   
      //D_tm GAMMA5, mu -> -mu
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spinout);
      dev_swapmu <<<1,1>>> ();
  
    //printf("Unbinding texture of spinorfield\n");
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    cudaMemcpy(spin4, spinout,size2, cudaMemcpyDeviceToDevice);
    //printf("Rebinding texture to spinorfield\n");
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
      
      //D_tm
      dev_tm_dirac_kappa<<<griddim3, blockdim3 >>>(gf, spin4, spin3, dev_nn);
    
    // r = b - Ax
    cublasSscal (24*VOLUME, -1.0, (float *)spin3, 1);
    cublasSaxpy (24*VOLUME, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    cublasScopy (24*VOLUME, (const float *)spin3, 1, (float *)spin0, 1);
    
    //dev_skalarmult_add_assign_spinor_field<<<griddim2, blockdim2 >>>(spinin, -1.0, spin3, spin0);
   }//recalculate residue

 }//MAIN LOOP cg	
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
 if(rescalekappa == 1){  //want D^-1 rescaled by 2*kappa
  
//multiply with D^dagger
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
      dev_gamma5 <<<griddim2, blockdim2 >>> (spin1,spin4);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spin4,1);
    #endif
      dev_tm_dirac_kappa <<<griddim3, blockdim3 >>> (gf, spin4, spin3, dev_nn);
      dev_gamma5 <<<griddim2, blockdim2 >>>(spin3,spin1);
      dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif


 //go over to non-kappa, Ddagger = g5 D g5
 dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spin1,1.0/(scaleparam*scaleparam), spinout);  
 
  // times operator == source ?? 
  //dev_tm_dirac_kappa<<<griddim3, blockdim3 >>>(gf, spin3, spinout, nn_grid);
  }
  else{
   dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  }
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
}





// this is the eo version of the device cg inner solver 
// we invert the hermitean Q_{-} Q_{+}
extern "C" void dev_cg_eo(dev_su3_2v * gf,dev_spinor* spinin, dev_spinor* spinout, 
dev_spinor* spin0, dev_spinor* spin1, dev_spinor* spin2, dev_spinor* spin3, dev_spinor* spin4, int *grid, int * nn_grid, REAL* output,REAL* erg, int xsize, int ysize, int rescalekappa, REAL epsfinal){
 
 
 REAL host_alpha, host_beta, host_dotprod, host_rk, sourcesquarenorm;
 REAL * dotprod, * dotprod2, * rk, * alpha, *beta;
 
 
 
 int i, gridsize;
 int maxit = max_innersolver_it;
 REAL eps = (REAL) innersolver_precision;
 int N_recalcres = 40; // after N_recalcres iterations calculate r = A x_k - b
 
 cudaError_t cudaerr;

 dim3 blockdim(1,1);
 dim3 blockdim2(128,1,1);
 if( VOLUME/2 >= 128){
   gridsize =VOLUME/2/128;
 }
 else{
   gridsize=1;
 }
 dim3 griddim2(gridsize,1,1);

 
 //if((VOLUME/2)%BLOCK != 0){
 //  printf("Error: VOLUME/2 is not a multiple of BLOCK. Aborting...\n");
 //  exit(100);
 //}
 int blockdim3=BLOCK;
 if( VOLUME/2 >= BLOCK){
   gridsize = (int)(VOLUME/2/BLOCK) + 1;
 }
 else{
   gridsize=1;
 }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
 
 // for dev_mul_one_pm...
 int blockdim4=BLOCK2;
 if( VOLUME/2 >= BLOCK2){
   gridsize = (int)(VOLUME/2/BLOCK2) + 1;
 }
 else{
   gridsize=1;
 }
 int griddim4=gridsize;  
 
 
 size_t size2 = sizeof(float4)*6*VOLUME/2;
 
 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(gf);
 #endif
 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  he_cg_init<<< 1, 1 >>> (grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3);
  // BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)
 
 // Init x,p,r for k=0
 // Allocate some numbers for host <-> device interaction
 cudaMalloc((void **) &dotprod, sizeof(REAL));
 cudaMalloc((void **) &dotprod2, sizeof(REAL));
 cudaMalloc((void **) &rk, sizeof(REAL));
 cudaMalloc((void **) &alpha, sizeof(REAL));
 cudaMalloc((void **) &beta, sizeof(REAL));
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 
 
 //init blas
 cublasInit();
 printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
 printf("have initialized cublas\n");
 
 
 // go over to kappa (if wanted)
 REAL scaleparam = sqrt(1.0/(2.0 * (REAL)hostkappa));
 printf("1/2kappa = %.8f\n",scaleparam);
 //dev_skalarmult_spinor_field<<<griddim2, blockdim2 >>>(spinin,scaleparam, spin1);
 //dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1, spinin);
 
 
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin0);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin1); // x_0 = 0
 dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spinin, spin2);
 dev_zero_spinor_field<<<griddim2, blockdim2 >>>(spin3);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
 
 

 //relative precision -> get initial residue
 sourcesquarenorm = cublasSdot (24*VOLUME/2, (const float *)spinin, 1, (const float *)spinin, 1);
 host_rk = sourcesquarenorm; //for use in main loop
 printf("Squarenorm Source:\t%.8e\n", sourcesquarenorm);
 printf("%s\n", cudaGetErrorString(cudaGetLastError()));
 
  printf("Entering cg-loop\n");
 for(i=0;i<maxit;i++){ //MAIN LOOP
  
  // Q_{-}Q{+}
  dev_Qtm_pm_psi(spin2, spin3, griddim3, blockdim3, griddim4, blockdim4);
  if((cudaerr=cudaGetLastError()) != cudaSuccess){
    printf("%s\n", cudaGetErrorString(cudaerr));
    exit(200);
  }
  
  
  
 //alpha
  host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin2, 1,
            (const float *) spin3, 1);
  host_alpha = (host_rk / host_dotprod); // alpha = r*r/ p M p
   
 //r(k+1)
 cublasSaxpy (24*VOLUME/2,-1.0*host_alpha, (const float *) spin3, 1, (float *) spin0, 1);  

 //x(k+1);
 cublasSaxpy (24*VOLUME/2, host_alpha, (const float *) spin2,  1, (float *) spin1, 1);

 printf("%s\n", cudaGetErrorString(cudaGetLastError()));

  //Abbruch?
  host_dotprod = cublasSdot (24*VOLUME/2, (const float *) spin0, 1,(const float *) spin0, 1);
  
 if (((host_dotprod <= eps*sourcesquarenorm) && (i > maxit / 2) ) || ( host_dotprod <= epsfinal/2.)){//error-limit erreicht (epsfinal/2 sollte ausreichen um auch in double precision zu bestehen)
   break; 
 }
  printf("iter %d: err = %.8e\n", i, host_dotprod);
  
 //beta
 host_beta =host_dotprod/host_rk;
 //p(k+1)
 cublasSscal (24*VOLUME/2, host_beta, (float *)spin2, 1);
 cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spin0,  1, (float *) spin2, 1);

 host_rk = host_dotprod;
 
 // recalculate residue frome r = b - Ax
 if(((i+1) % N_recalcres) == 0){
    // r_(k+1) = Ax -b 
    printf("Recalculating residue\n");
    
    // D Ddagger   --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!
      
    // Q_{-}Q{+}
    dev_Qtm_pm_psi(spin1, spin3, griddim3, blockdim3, griddim4, blockdim4);
      
        
    
    // r = b - Ax
    cublasSscal (24*VOLUME/2, -1.0, (float *)spin3, 1);
    cublasSaxpy (24*VOLUME/2, 1.0, (const float *) spinin,  1, (float *) spin3, 1);
    cublasScopy (24*VOLUME/2, (const float *)spin3, 1, (float *)spin0, 1);
    //dev_skalarmult_add_assign_spinor_field<<<griddim2, blockdim2 >>>(spinin, -1.0, spin3, spin0);
   }//recalculate residue

 }//MAIN LOOP cg	
  
  
  printf("Final residue: %.6e\n",host_dotprod);
  // x_result = spin1 !
  
  //no multiplication with D^{dagger} here and no return to non-kappa basis as in dev_cg!
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(spin1,spinout);
  
  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(dotprod);
  cudaFree(dotprod2);
  cudaFree(rk);
  cudaFree(alpha);
  cudaFree(beta);
  cublasShutdown();
}












//initialize nearest-neighbour table for gpu
void initnn(){
  int t,x,y,z,pos;
  for(t=0;t<T;t++){
   for(x=0; x<LX; x++){
    for(y=0; y<LY; y++){
     for(z=0; z<LZ; z++){   
          pos= z + LZ*(y + LY*(x + LX*t));
          //plus direction
          nn[8*pos+0] = z + LZ*(y + LY*(x + LX*((t+1)%T)));
          nn[8*pos+1] = z + LZ*(y + LY*((x+1)%LX + LX*t));
          nn[8*pos+2] = z + LZ*((y+1)%LY + LY*(x + LX*t));
          nn[8*pos+3] = (z+1)%LZ + LX*(y + LY*(x + LX*t));
          //minus direction
          if(t==0){
            nn[8*pos+4] = z + LZ*(y + LY*(x + LX*((T-1))));
          }
          else{
            nn[8*pos+4] = z + LZ*(y + LY*(x + LX*((t-1))));
          }
          if(x==0){
            nn[8*pos+5] = z + LZ*(y + LY*((LX-1) + LX*t));
          }
          else{
            nn[8*pos+5] = z + LZ*(y + LY*((x-1) + LX*t));
          }
          if(y==0){
            nn[8*pos+6] = z + LZ*((LY-1) + LY*(x + LX*t));
          }
          else{
            nn[8*pos+6] = z + LZ*((y-1) + LY*(x + LX*t));
          }
          if(z==0){
            nn[8*pos+7] = (LZ-1) + LZ*(y + LY*(x + LX*t));
          }
          else{
            nn[8*pos+7] = (z-1) + LZ*(y + LY*(x + LX*t));
          }          
        }
      }
    } 
  }
}





//initialize nearest-neighbour table for gpu with even-odd enabled
//init_nn must have been called before for initialization of nn
void initnn_eo(){
  int x,y,z,t,ind,nnpos,j;
  int evenpos=0;
  int oddpos=0;
  for(t=0;t<T;t++){
    for(x=0;x<LX;x++){
      for(y=0;y<LY;y++){
        for(z=0;z<LZ;z++){
          ind = g_ipt[t][x][y][z];
          
          if(((t+x+y+z)%2 == 0)){
            nnpos = g_lexic2eosub[ind];
            for(j=0;j<4;j++){
              nn_eo[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];
            }
            for(j=0;j<4;j++){
              nn_eo[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
            eoidx_even[evenpos] = ind;
            evenpos++;
          }
          else{
            nnpos = g_lexic2eosub[ind];
            for(j=0;j<4;j++){
              nn_oe[8*nnpos+j] = g_lexic2eosub[ g_iup[ind][j] ];
            }
            for(j=0;j<4;j++){
              nn_oe[8*nnpos+4+j] = g_lexic2eosub[ g_idn[ind][j] ];
            }
            eoidx_odd[oddpos] = ind;
            oddpos++;
          }
        }
      }
    }
  }
}




// show the nn table eo
void shownn_eo(){
  int i,pos;
  printf("eo part\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
       for(i=0;i<8;i++){
          printf("%d ",nn_eo[8*pos+i]);
          //lptovec(nn[8*pos+i]);
        }
        printf("\n");
    }
  printf("oe part\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
       for(i=0;i<8;i++){
          printf("%d ",nn_oe[8*pos+i]);
          //lptovec(nn[8*pos+i]);
        }
        printf("\n");
    }
    
  printf("site index even\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
          printf("%d ",eoidx_even[pos]);
          //lptovec(nn[8*pos+i]);
        printf("\n");
  }

  printf("site index odd\n");
  for(pos=0;pos<VOLUME/2;pos++){ 
       printf("p=%d\t", pos);
          printf("%d ",eoidx_odd[pos]);
          //lptovec(nn[8*pos+i]);
        printf("\n");
  }
  printf("checking forward even\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_oe[8*nn_eo[8*pos+i]+4+i]);
    }
  }

  printf("checking backward even\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_oe[8*nn_eo[8*pos+4+i]+i]);
    }
  }

  printf("checking forward odd\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_eo[8*nn_oe[8*pos+i]+4+i]);
    }
  }

  printf("checking backward odd\n");
  for(pos=0;pos<VOLUME/2;pos++){
    for(i=0;i<4;i++){
      printf("%d = %d\n",pos, nn_eo[8*nn_oe[8*pos+4+i]+i]);
    }
  }
}


void show_su3(su3 gf1){
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c00.re,
   					gf1.c00.im,
   					gf1.c01.re,
   					gf1.c01.im,
   					gf1.c02.re,
   					gf1.c02.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c10.re,
   					gf1.c10.im,
   					gf1.c11.re,
   					gf1.c11.im,
   					gf1.c12.re,
   					gf1.c12.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1.c20.re,
   					gf1.c20.im,
   					gf1.c21.re,
   					gf1.c21.im,
   					gf1.c22.re,
   					gf1.c22.im
   ); 
}


void show_dev_su3(dev_su3 gf1){
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[0][0].re,
   					gf1[0][0].im,
   					gf1[0][1].re,
   					gf1[0][1].im,
   					gf1[0][2].re,
   					gf1[0][2].im
   );   
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[1][0].re,
   					gf1[1][0].im,
   					gf1[1][1].re,
   					gf1[1][1].im,
   					gf1[1][2].re,
   					gf1[1][2].im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[2][0].re,
   					gf1[2][0].im,
   					gf1[2][1].re,
   					gf1[2][1].im,
   					gf1[2][2].re,
   					gf1[2][2].im
   ); 

}


void lptovec(int k){
  int L3 = L*L*L;
  int L2 = L*L;
  int x0,x1,x2,x3;
  x0 = k/L3;
  k = k-x0*L3; 
  x3 = k/L2;
  k = k-x3*L2;
  x2 = k/L;
  k = k-x2*L;
  x1 = k;
  printf("%d,%d,%d,%d;  ",x0,x3,x2,x1);
}


// show nn table 
void shownn(){
  int t,x,y,z,i,pos;
  int lx,ly,lz,lt;
    lx = LX;
    ly = LY;
    lz = LZ;
    lt =T;  
  for(t=0;t<lt;t++){ 
    for(x=0; x<lx; x++){
      for(y=0; y<ly; y++){
        for(z=0; z<lz; z++){
          pos= z + lz*(y + ly*(x + lx*t));
          printf("p=%d\t", pos);
          for(i=0;i<8;i++){
            printf("%d ",nn[8*pos+i]);
            //lptovec(nn[8*pos+i]);
          }
          printf("\n");
          //compare with geometry fields of hmc
          //might NOT WORK for even-odd? What are geometry indices in case of even-odd?
          printf("%d: %d %d %d %d %d %d %d %d\n",g_ipt[t][x][y][z],g_iup[pos][0],g_iup[pos][1],g_iup[pos][2],g_iup[pos][3],g_idn[pos][0],g_idn[pos][1],g_idn[pos][2],g_idn[pos][3]);
        }
      }
    }
  }
}




// get 2 first rows of gf float4 type
//  
//
void su3to2vf4(su3** gf, dev_su3_2v* h2d_gf){
  int i,j;
  for (i=0;i<VOLUME;i++){
   for(j=0;j<4;j++){
   //first row
    h2d_gf[3*(4*i+j)].x = (REAL) gf[i][j].c00.re;
    h2d_gf[3*(4*i+j)].y = (REAL) gf[i][j].c00.im;
    h2d_gf[3*(4*i+j)].z = (REAL) gf[i][j].c01.re;
    h2d_gf[3*(4*i+j)].w = (REAL) gf[i][j].c01.im;
    h2d_gf[3*(4*i+j)+1].x = (REAL) gf[i][j].c02.re;
    h2d_gf[3*(4*i+j)+1].y = (REAL) gf[i][j].c02.im;      
   //second row
    h2d_gf[3*(4*i+j)+1].z = (REAL) gf[i][j].c10.re;
    h2d_gf[3*(4*i+j)+1].w = (REAL) gf[i][j].c10.im;
    h2d_gf[3*(4*i+j)+2].x = (REAL) gf[i][j].c11.re;
    h2d_gf[3*(4*i+j)+2].y = (REAL) gf[i][j].c11.im;
    h2d_gf[3*(4*i+j)+2].z = (REAL) gf[i][j].c12.re;
    h2d_gf[3*(4*i+j)+2].w = (REAL) gf[i][j].c12.im;      
  } 
 }
}




// bring gf into the form
// a2 a3, theta_a1, theta_c1, b1
// 
void su3to8(su3** gf, dev_su3_8* h2d_gf){
  int i,j;
  for (i=0;i<VOLUME;i++){
   for(j=0;j<4;j++){
   // a2, a3
    h2d_gf[2*(4*i+j)].x = (REAL) gf[i][j].c01.re;
    h2d_gf[2*(4*i+j)].y = (REAL) gf[i][j].c01.im;
    h2d_gf[2*(4*i+j)].z = (REAL) gf[i][j].c02.re;
    h2d_gf[2*(4*i+j)].w = (REAL) gf[i][j].c02.im;
    
   // theta_a1, theta_c1
   // use atan2 for this: following the reference, atan2 should give an angle -pi < phi < +pi  
   h2d_gf[2*(4*i+j)+1].x = (REAL)( atan2((REAL) gf[i][j].c00.im,(REAL) gf[i][j].c00.re ));
   h2d_gf[2*(4*i+j)+1].y = (REAL) ( atan2((REAL) gf[i][j].c20.im,(REAL)gf[i][j].c20.re ));
     
   // b1
    h2d_gf[2*(4*i+j)+1].z = (REAL) gf[i][j].c10.re ;
    h2d_gf[2*(4*i+j)+1].w = (REAL) gf[i][j].c10.im ;
     
  } 
 }
}






// this is to reconstruct the gf on the host from 2 rows of the link
// may be used for tests
void reconstructgf_2v (dev_su3* gf){
  complex help1;
  complex help2;
  //third row from cconj(cross product of first and second row)
  _mult_assign_complex(help1,(*gf)[0][1],(*gf)[1][2]);
  _mult_assign_complex(help2,(*gf)[0][2],(*gf)[1][1]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][0].re = help1.re;
  (*gf)[2][0].im = help1.im;
  
  _mult_assign_complex(help1,(*gf)[0][2],(*gf)[1][0]);
  _mult_assign_complex(help2,(*gf)[0][0],(*gf)[1][2]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][1].re = help1.re;
  (*gf)[2][1].im = help1.im;
  
  _mult_assign_complex(help1,(*gf)[0][0],(*gf)[1][1]);
  _mult_assign_complex(help2,(*gf)[0][1],(*gf)[1][0]);
  _diff_complex(help1,help2);
  help1.im = -help1.im;
  (*gf)[2][2].re = help1.re;
  (*gf)[2][2].im = help1.im;
  return;
}







__global__ void dev_check_gauge_reconstruction_8(dev_su3_2v* gf, int pos, dev_su3 * outgf1, dev_su3* outgf2){
  dev_reconstructgf_8texref (gf,pos, outgf1);
  dev_reconstructgf_8texref_dagger (gf,pos, outgf2);
}







void check_gauge_reconstruction_8(su3 ** gf1, dev_su3_2v * gf2, int ind1, int mu){
  dev_su3 * reconst_g , * reconst_g_dagger;
  dev_su3  result, result_dagger;
   printf("Checking 8 paramater reconstruction of gauge field:\n");  
  su3 gfdagger;
  #ifdef USETEXTURE
    bind_texture_gf(gf2);
  #endif
  printf("\n");
  size_t cpsize = sizeof(dev_su3); // parallel in t and z direction
  cudaMalloc((void **) &reconst_g, cpsize); 
  cudaMalloc((void **) &reconst_g_dagger, cpsize); 
  
  show_su3(gf1[ind1][mu]);
  printf("\n");
  
  dev_check_gauge_reconstruction_8  <<< 1 , 1 >>> (dev_gf,4*ind1 + mu, reconst_g, reconst_g_dagger);
  cudaMemcpy(&result, reconst_g, cpsize, cudaMemcpyDeviceToHost);
  cudaMemcpy(&result_dagger, reconst_g_dagger, cpsize, cudaMemcpyDeviceToHost);

  show_dev_su3(result);
  printf("\n");
  
  _su3_dagger(gfdagger,gf1[ind1][mu]);
  show_su3(gfdagger);
  printf("\n");
  show_dev_su3(result_dagger);



  #ifdef USETEXTURE
   unbind_texture_gf();
  #endif
  cudaFree(reconst_g);
}





// this is to reconstruct the gf on the host from 2 rows of the link
// may be used for tests
void reconstructgf_8 (dev_su3_8 * h2d_gf, dev_su3* gf){

  float4 gfin;
  REAL N, one_over_N, help;
  complex p1,p2, chelp1, chelp2, chelp3, chelpconj, chelpconj2;
  
  gfin = h2d_gf[0];
  // read a2 a3
  (*gf)[0][1].re = gfin.x;
  (*gf)[0][1].im = gfin.y;
  (*gf)[0][2].re = gfin.z;
  (*gf)[0][2].im = gfin.w;  
 
  help = gfin.x*gfin.x + gfin.y*gfin.y + gfin.z*gfin.z + gfin.w*gfin.w; // use later on
  N = sqrt(help);
  one_over_N = 1.0f/N;
  
  // read theta_a1, theta_c1, b1
  gfin = h2d_gf[1];
  
  // reconstruct a1
  help = sqrt(1.0f - help);
  (*gf)[0][0].re = help*cos(gfin.x);
  (*gf)[0][0].im = help*sin(gfin.x);
  
  // assign b1
  (*gf)[1][0].re = gfin.z;
  (*gf)[1][0].im = gfin.w;
  
  // p2 = 1/N b1
  p2.re = one_over_N*(*gf)[1][0].re;
  p2.im = one_over_N*(*gf)[1][0].im;  


  // reconstruct c1
  help = sqrt(1.0f - 
              (*gf)[0][0].re * (*gf)[0][0].re - (*gf)[0][0].im * (*gf)[0][0].im - 
              (*gf)[1][0].re * (*gf)[1][0].re - (*gf)[1][0].im * (*gf)[1][0].im
          );
  (*gf)[2][0].re = help*cos(gfin.y);
  (*gf)[2][0].im = help*sin(gfin.y);

  
  // p1 = 1/N*cconj(c1)
  p1.re = one_over_N*(*gf)[2][0].re;
  p1.im = - one_over_N*(*gf)[2][0].im;
  
  
  float temp = p1.re*p1.re + p1.im*p1.im + p2.re*p2.re + p2.im*p2.im;
  printf("p1**2 + p2**2 = %f\n", temp);
  
  
  // calculate b2
  _complex_conj(chelpconj, (*gf)[0][2] );
  _mult_assign_complex(chelp1, p1, chelpconj   );
  _complex_conj(chelpconj, (*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj , (*gf)[0][1] ); 
  _mult_assign_complex(chelp2, p2, chelp3);
  _add_complex(chelp1, chelp2);
  _mult_real((*gf)[1][1], chelp1, -one_over_N);

  
  // calculate b3
  _complex_conj(chelpconj, (*gf)[0][1] );
  _mult_assign_complex(chelp1, p1,  chelpconj  );
  _complex_conj(chelpconj, (*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj  , (*gf)[0][2] );  
  _mult_assign_complex(chelp2, p2, chelp3 );
  _diff_complex(chelp1, chelp2);
  _mult_real((*gf)[1][2],chelp1, one_over_N);

  
  // calculate c2
  _complex_conj(chelpconj, p2);
  _complex_conj(chelpconj2, (*gf)[0][2]);
  _mult_assign_complex(chelp1, chelpconj , chelpconj2 );
  _complex_conj(chelpconj,(*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj  , (*gf)[0][1] );
  _complex_conj(chelpconj2,p1);
  _mult_assign_complex(chelp2, chelpconj2  , chelp3);
  _diff_complex(chelp1, chelp2);
  _mult_real((*gf)[2][1],chelp1, one_over_N);
  
  
  // calculate c3
  _complex_conj(chelpconj, p2);
  _complex_conj(chelpconj2, (*gf)[0][1] );
  _mult_assign_complex(chelp1, chelpconj  , chelpconj2   );
  _complex_conj(chelpconj,(*gf)[0][0]);
  _mult_assign_complex(chelp3, chelpconj ,(*gf)[0][2]);
  _complex_conj(chelpconj,p1);
  _mult_assign_complex( chelp2, chelpconj  , chelp3 );
  _add_complex(chelp1, chelp2);
  _mult_real((*gf)[2][2], chelp1, -one_over_N);
                  
}




void showcompare_gf(int t, int x, int y, int z, int mu){
   int ind1 = g_ipt[t][x][y][z];
   su3 ** gf1 = g_gauge_field;
   
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c00.re,
   					gf1[ind1][mu].c00.im,
   					gf1[ind1][mu].c01.re,
   					gf1[ind1][mu].c01.im,
   					gf1[ind1][mu].c02.re,
   					gf1[ind1][mu].c02.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c10.re,
   					gf1[ind1][mu].c10.im,
   					gf1[ind1][mu].c11.re,
   					gf1[ind1][mu].c11.im,
   					gf1[ind1][mu].c12.re,
   					gf1[ind1][mu].c12.im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",gf1[ind1][mu].c20.re,
   					gf1[ind1][mu].c20.im,
   					gf1[ind1][mu].c21.re,
   					gf1[ind1][mu].c21.im,
   					gf1[ind1][mu].c22.re,
   					gf1[ind1][mu].c22.im
   );
   printf("\n\n");

   int ind2 =  z + LZ*(y + LY*(x + LX*t));
#ifdef GF_8
   printf("8-field:\t(%f,%f,%f,%f) (%f,%f,%f,%f)\n",
     h2d_gf[2*(4*ind2+mu)].x,
     h2d_gf[2*(4*ind2+mu)].y,
     h2d_gf[2*(4*ind2+mu)].z,
     h2d_gf[2*(4*ind2+mu)].w,
     h2d_gf[2*(4*ind2+mu)+1].x,
     h2d_gf[2*(4*ind2+mu)+1].y,
     h2d_gf[2*(4*ind2+mu)+1].z,
     h2d_gf[2*(4*ind2+mu)+1].w
   );
   dev_su3 help; 
   reconstructgf_8( &(h2d_gf[2*(4*ind2+mu)]) , &help );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[0][0].re,
   					help[0][0].im,
   					help[0][1].re,
   					help[0][1].im,
   					help[0][2].re,
   					help[0][2].im
   );   
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[1][0].re,
   					help[1][0].im,
   					help[1][1].re,
   					help[1][1].im,
   					help[1][2].re,
   					help[1][2].im
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[2][0].re,
   					help[2][0].im,
   					help[2][1].re,
   					help[2][1].im,
   					help[2][2].re,
   					help[2][2].im
   );   
   
#else
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",h2d_gf[3*(4*ind2+mu)].x,
   					h2d_gf[3*(4*ind2+mu)].y,
   					h2d_gf[3*(4*ind2+mu)].z,
   					h2d_gf[3*(4*ind2+mu)].w,
   					h2d_gf[3*(4*ind2+mu)+1].x,
   					h2d_gf[3*(4*ind2+mu)+1].y
   );
   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",h2d_gf[3*(4*ind2+mu)+1].z,
   					h2d_gf[3*(4*ind2+mu)+1].w,
   					h2d_gf[3*(4*ind2+mu)+2].x,
   					h2d_gf[3*(4*ind2+mu)+2].y,
   					h2d_gf[3*(4*ind2+mu)+2].z,
   					h2d_gf[3*(4*ind2+mu)+2].w
   );   
   
   dev_su3 help;
   
   help[0][0].re = h2d_gf[3*(4*ind2+mu)].x;
   help[0][0].im = h2d_gf[3*(4*ind2+mu)].y;
   help[0][1].re = h2d_gf[3*(4*ind2+mu)].z;
   help[0][1].im = h2d_gf[3*(4*ind2+mu)].w;

   help[0][2].re = h2d_gf[3*(4*ind2+mu)+1].x;
   help[0][2].im = h2d_gf[3*(4*ind2+mu)+1].y;
   help[1][0].re = h2d_gf[3*(4*ind2+mu)+1].z;
   help[1][0].im = h2d_gf[3*(4*ind2+mu)+1].w;
   
   help[1][1].re = h2d_gf[3*(4*ind2+mu)+2].x;
   help[1][1].im = h2d_gf[3*(4*ind2+mu)+2].y;
   help[1][2].re = h2d_gf[3*(4*ind2+mu)+2].z;
   help[1][2].im = h2d_gf[3*(4*ind2+mu)+2].w;   
   
   reconstructgf_2v (&help); 

   printf("(%f,%f)\t(%f,%f)\t(%f,%f)\n",help[2][0].re,
   					help[2][0].im,
   					help[2][1].re,
   					help[2][1].im,
   					help[2][2].re,
   					help[2][2].im
   );
#endif 
}








// convert spinor to double 
void convert2double_spin(dev_spinor* spin, spinor* h2d){
  int i,Vol;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
  
        h2d[i].s0.c0.re = (double) spin[6*i+0].x;
        h2d[i].s0.c0.im = (double) spin[6*i+0].y;
        h2d[i].s0.c1.re = (double) spin[6*i+0].z;
        h2d[i].s0.c1.im = (double) spin[6*i+0].w;
        
        h2d[i].s0.c2.re = (double) spin[6*i+1].x;
        h2d[i].s0.c2.im = (double) spin[6*i+1].y;
        h2d[i].s1.c0.re = (double) spin[6*i+1].z;
        h2d[i].s1.c0.im = (double) spin[6*i+1].w;   
        
        h2d[i].s1.c1.re = (double) spin[6*i+2].x;
        h2d[i].s1.c1.im = (double) spin[6*i+2].y;
        h2d[i].s1.c2.re = (double) spin[6*i+2].z;
        h2d[i].s1.c2.im = (double) spin[6*i+2].w;  
        
        h2d[i].s2.c0.re = (double) spin[6*i+3].x;
        h2d[i].s2.c0.im = (double) spin[6*i+3].y;
        h2d[i].s2.c1.re = (double) spin[6*i+3].z;
        h2d[i].s2.c1.im = (double) spin[6*i+3].w;  
        
        h2d[i].s2.c2.re = (double) spin[6*i+4].x;
        h2d[i].s2.c2.im = (double) spin[6*i+4].y;
        h2d[i].s3.c0.re = (double) spin[6*i+4].z;
        h2d[i].s3.c0.im = (double) spin[6*i+4].w; 
        
        h2d[i].s3.c1.re = (double) spin[6*i+5].x;
        h2d[i].s3.c1.im = (double) spin[6*i+5].y;
        h2d[i].s3.c2.re = (double) spin[6*i+5].z;
        h2d[i].s3.c2.im = (double) spin[6*i+5].w; 
        
  }
}





// convert spinor to REAL4 (float4, double4) 
void convert2REAL4_spin(spinor* spin, dev_spinor* h2d){
  int i,Vol;
  if(even_odd_flag){
    Vol = VOLUME/2;
  }
  else{
    Vol = VOLUME;
  }
  for (i=0;i<Vol;i++){
    
        h2d[6*i+0].x = (REAL) spin[i].s0.c0.re;
        h2d[6*i+0].y = (REAL) spin[i].s0.c0.im;
        h2d[6*i+0].z = (REAL) spin[i].s0.c1.re;
        h2d[6*i+0].w = (REAL) spin[i].s0.c1.im;
        
        h2d[6*i+1].x = (REAL) spin[i].s0.c2.re;
        h2d[6*i+1].y = (REAL) spin[i].s0.c2.im;
        h2d[6*i+1].z = (REAL) spin[i].s1.c0.re;
        h2d[6*i+1].w = (REAL) spin[i].s1.c0.im;
        
        h2d[6*i+2].x = (REAL) spin[i].s1.c1.re;
        h2d[6*i+2].y = (REAL) spin[i].s1.c1.im;
        h2d[6*i+2].z = (REAL) spin[i].s1.c2.re;
        h2d[6*i+2].w = (REAL) spin[i].s1.c2.im;
        
        h2d[6*i+3].x = (REAL) spin[i].s2.c0.re;
        h2d[6*i+3].y = (REAL) spin[i].s2.c0.im;
        h2d[6*i+3].z = (REAL) spin[i].s2.c1.re;
        h2d[6*i+3].w = (REAL) spin[i].s2.c1.im;
        
        h2d[6*i+4].x = (REAL) spin[i].s2.c2.re;
        h2d[6*i+4].y = (REAL) spin[i].s2.c2.im;
        h2d[6*i+4].z = (REAL) spin[i].s3.c0.re;
        h2d[6*i+4].w = (REAL) spin[i].s3.c0.im;
        
        h2d[6*i+5].x = (REAL) spin[i].s3.c1.re;
        h2d[6*i+5].y = (REAL) spin[i].s3.c1.im;
        h2d[6*i+5].z = (REAL) spin[i].s3.c2.re;
        h2d[6*i+5].w = (REAL) spin[i].s3.c2.im;
    
  }
}





void init_mixedsolve(su3** gf){
cudaError_t cudaerr;

   // get number of devices
   if(havedevice == 0){
     int ndev = find_devices();
	   if(ndev == 0){
	       fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
	       exit(300);
	    }
	 // try to set active device to device_num given in input file
	    if(device_num < ndev){
	     printf("Setting active device to: %d\n", device_num);
	     cudaSetDevice(device_num);
	    }
	    else{
	      fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
	      exit(301);
	    }
	    if((cudaerr=cudaGetLastError())!=cudaSuccess){
	    printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
	    exit(302);
	    }
    havedevice = 1;
    }
  #ifdef GF_8
  /* allocate 8 floats of gf = 2*4*VOLUME float4's*/
  printf("Using GF 8 reconstruction\n");
  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8);
  #else
  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  printf("Using GF 12 reconstruction\n");
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 
  #endif
  
  #ifdef USETEXTURE
    printf("Using texture references\n");
  #else
    printf("NOT using texture references\n");
  #endif
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  #ifdef GF_8
  h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to8(gf,h2d_gf);  
  #else
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4(gf,h2d_gf);
  #endif
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);


//grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn = (int *) malloc(nnsize);
  cudaMalloc((void **) &dev_nn, nnsize);
  
  initnn();
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
  
  //free again
  free(nn);


// Spinors
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); /* float4 */

  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  
  cudaMalloc((void **) &dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &dev_spin4, dev_spinsize);
  cudaMalloc((void **) &dev_spin5, dev_spinsize);
  cudaMalloc((void **) &dev_spinin, dev_spinsize);
  cudaMalloc((void **) &dev_spinout, dev_spinsize);
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  
  
  output_size = LZ*T*sizeof(float); // parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);   // output array
  float * host_output = (float*) malloc(output_size);

  int grid[5];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME;
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
}






void init_mixedsolve_eo(su3** gf){
cudaError_t cudaerr;
  dev_complex help;

  if(havedevice == 0){
   // get number of devices
     int ndev = find_devices();
	   if(ndev == 0){
	       fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
	       exit(300);
	    }
	 // try to set active device to device_num given in input file
	    if(device_num < ndev){
	     printf("Setting active device to: %d\n", device_num);
	     cudaSetDevice(device_num);
	    }
	    else{
	      fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
	      exit(301);
	    }
	    if((cudaerr=cudaGetLastError())!=cudaSuccess){
	    printf("Error in init_mixedsolve_eo(): Could not set active device. Aborting...\n");
	    exit(302);
	    }  
   havedevice=1;
  }
  #ifdef GF_8
  /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
  printf("Using GF 8 reconstruction\n");
  size_t dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8); 
  #else
  /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
  printf("Using GF 12 reconstruction\n");
  size_t dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 
  #endif
  
  #ifdef USETEXTURE
    printf("Using texture references\n");
  #else
    printf("NOT using texture references\n");
  #endif
  
  if((cudaerr=cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of gauge field failed. Aborting...\n");
    exit(200);
  }   // Allocate array on device
  else{
    printf("Allocated gauge field on device\n");
  }  
  
  
  
  #ifdef GF_8
  h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to8(gf,h2d_gf);
  #else
  h2d_gf = (dev_su3_2v *)malloc(dev_gfsize); // Allocate REAL conversion gf on host
  su3to2vf4(gf,h2d_gf);
  #endif
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);



//grid 
  size_t nnsize = 8*VOLUME*sizeof(int);
  nn = (int *) malloc(nnsize);
  
  //nn grid for even-odd
  nn_eo = (int *) malloc(nnsize/2);
  nn_oe = (int *) malloc(nnsize/2);
  
  cudaMalloc((void **) &dev_nn, nnsize);
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);
  
  
  size_t idxsize = VOLUME/2*sizeof(int);
  eoidx_even = (int *) malloc(idxsize);
  eoidx_odd = (int *) malloc(idxsize);
  cudaMalloc((void **) &dev_eoidx_even, idxsize);
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);
  
  initnn();
  initnn_eo();
  //shownn_eo();
  
  //shownn();
  //showcompare_gf(T-1, LX-1, LY-1, LZ-1, 3);
  //check_gauge_reconstruction_8(gf, dev_gf, 0, 0);
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);
  
  //free again
  free(eoidx_odd);
  free(eoidx_even);
  free(nn_oe);
  free(nn_eo);
  free(nn);
  
// Spinors
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); /* float4 */

  if((void*)(h2d_spin = (dev_spinor *)malloc(dev_spinsize)) == NULL){
    printf("Could not allocate memory for h2d_spin. Aborting...\n");
    exit(200);
  } // Allocate float conversion spinor on host
  
  cudaMalloc((void **) &dev_spin1, dev_spinsize);   // Allocate array spin1 on device
  cudaMalloc((void **) &dev_spin2, dev_spinsize);   // Allocate array spin2 on device
  cudaMalloc((void **) &dev_spin3, dev_spinsize);   // Allocate array spin3 on device
  cudaMalloc((void **) &dev_spin4, dev_spinsize);
  cudaMalloc((void **) &dev_spin5, dev_spinsize);
  cudaMalloc((void **) &dev_spinin, dev_spinsize);
  cudaMalloc((void **) &dev_spinout, dev_spinsize);
  
  cudaMalloc((void **) &dev_spin_eo1, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo2, dev_spinsize);
  
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve(): Memory allocation of spinor fields failed. Aborting...\n");
    exit(200);
  }
  else{
    printf("Allocated spinor fields on device\n");
  }
  
  
  output_size = LZ*T*sizeof(float); // parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);   // output array
  float * host_output = (float*) malloc(output_size);

  int grid[5];
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2; 
  // dev_VOLUME is half of VOLUME for eo
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
}



void finalize_mixedsolve(){

  cudaFree(dev_spin1);
  cudaFree(dev_spin2);
  cudaFree(dev_spin3);
  cudaFree(dev_spin4);
  cudaFree(dev_spin5);
  cudaFree(dev_spinin);
  cudaFree(dev_spinout);
  cudaFree(dev_gf);
  cudaFree(dev_grid);
  cudaFree(dev_output);
  cudaFree(dev_nn);
  
  if(even_odd_flag){
    cudaFree(dev_spin_eo1);
    cudaFree(dev_spin_eo2);
    cudaFree(dev_eoidx_even);
    cudaFree(dev_eoidx_odd);
    cudaFree(dev_nn_eo);
    cudaFree(dev_nn_oe);
  }
  
  
  
  free(h2d_spin);
  free(h2d_gf);
}







extern "C" int mixed_solve (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec,const int N){
  
  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter;
  
  size_t dev_spinsize = 6*VOLUME * sizeof(dev_spinor); // float4 
  init_mixedsolve(g_gauge_field);
  
  // Start timer
  assert((start = clock())!=-1);
  
  rk = square_norm(Q, N, 1);
  sourcesquarenorm = rk; // for relative precision
  assign(g_spinor_field[DUM_SOLVER],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(g_spinor_field[DUM_SOLVER+1],  N);//spin2 = x_k
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],  N);
  printf("The VOLUME is: %d\n",N);
  
  
  
for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision Dirac-Op...\n");
   
   Q_pm_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
   diff(g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER+3],N);
    // r_k = b - D x_k
   
   rk = square_norm(g_spinor_field[DUM_SOLVER], N, 0);
  
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
    finalize_mixedsolve();
    return(iter*max_innersolver_it);  // MAYBE ONE SHOULD KEEP TRACK OF REAL INNER SOLVER STEPS
   }
   

  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   dev_cg(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0);
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
   
  
   // copy back
   cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
   convert2double_spin(h2d_spin, g_spinor_field[DUM_SOLVER+2]);
   
   add(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+2],N);
   // x_(k+1) = x_k + p_k
   
   outercount ++;
    
}// outer loop 

     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with D^dagger
     Q_minus_psi_gpu(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
  

    stop = clock();
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  return(-1);
}




void dummy (dev_spinor* a, dev_spinor* b){

}


void benchmark(spinor * const Q){
  
  double timeelapsed = 0.0;
  clock_t start, stop;
  int i;
  
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
  convert2REAL4_spin(Q,h2d_spin);
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  assert((start = clock())!=-1);

 #ifdef USETEXTURE
  //Bind texture gf
  bind_texture_gf(dev_gf);
 #endif

 //Initialize some stuff
  printf("mu = %f\n", g_mu);
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex)) ;
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)) ; 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex)) ;  
  
  
  int blockdim3=BLOCK;
  int gridsize;
  if( VOLUME/2 >= BLOCK){
    gridsize = (int)(VOLUME/2/BLOCK) + 1;
  }
  else{
    gridsize=1;
  }
 printf("gridsize = %d\n", gridsize);
 int griddim3=gridsize; 
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0,h1,h2,h3); 
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));
  printf("Applying H 1000 times\n");
  for(i=0; i<1000; i++){
      #ifdef USETEXTURE
       bind_texture_spin(dev_spinin,1);
      #endif
       //bind_texture_nn(dev_nn_eo);
      //cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
      dev_Hopping_Matrix<<<griddim3, blockdim3>>>
             (dev_gf, dev_spinin, dev_spin_eo1, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0); //dev_spin_eo1 == even -> 0
       //unbind_texture_nn();
    #ifdef USETEXTURE             
     unbind_texture_spin(1);
    #endif

    #ifdef USETEXTURE
     bind_texture_spin(dev_spin_eo1,1);
    #endif
  //bind_texture_nn(dev_nn_oe);
   // cudaFuncSetCacheConfig(dev_Hopping_Matrix, cudaFuncCachePreferL1);
    dev_Hopping_Matrix<<<griddim3, blockdim3>>>
            (dev_gf, dev_spin_eo1, dev_spinin, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1); 
  //unbind_texture_nn();
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif

  }  
  printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
  printf("Done\n"); 
  
  assert((stop = clock())!=-1);
  timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
  // x2 because 2x Hopping per iteration
  double benchres = 1400.0*2*(VOLUME/2)* 1000 / timeelapsed / 1.0e9;
  printf("Benchmark: %f Gflops\n", benchres); 
   
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
}





extern "C" int mixed_solve_eo (spinor * const P, spinor * const Q, const int max_iter, 
	   double eps, const int rel_prec, const int N){

  // source in Q, initial solution in P (not yet implemented)
  double rk;
  int outercount=0;
  clock_t start, stop, startinner, stopinner; 
  double timeelapsed = 0.0;
  double sourcesquarenorm;
  int iter, retval;
  

  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor); // float4 even-odd !
    
  init_mixedsolve_eo(g_gauge_field);
  
  /*
  // small benchmark
    assign(g_spinor_field[DUM_SOLVER],Q,N);
    benchmark(g_spinor_field[DUM_SOLVER]);
  // end small benchmark
  
  exit(100);
  */
 
 


  // Start timer
  assert((start = clock())!=-1);
  rk = square_norm(Q, N, 1);
  sourcesquarenorm=rk; // for relative prec
  double finaleps;
  if(rel_prec == 1){
    finaleps = eps * sourcesquarenorm;
  }
  else{
    finaleps = eps;
  }
  assign(g_spinor_field[DUM_SOLVER],Q,N);
  printf("Initial residue: %.16e\n",rk);
  zero_spinor_field(g_spinor_field[DUM_SOLVER+1],  N);//spin2 = x_k
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],  N);
  printf("The VOLUME/2 is: %d\n",N);
  
for(iter=0; iter<max_iter; iter++){

   printf("Applying double precision EO Dirac-Op Q_{-}Q{+}...\n");
   
   Qtm_pm_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+2]);
   diff(g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER],g_spinor_field[DUM_SOLVER+3],N);
    // r_k = b - D x_k
   
   rk = square_norm(g_spinor_field[DUM_SOLVER], N, 0);
   #ifdef GF_8
    if(isnan(rk)){
      fprintf(stderr, "Error in mixed_solve_eo: Residue is NaN.\n  May happen with GF 8 reconstruction. Aborting ...\n");
      exit(200);
    }
   #endif
   
   printf("Residue after %d inner solver iterations: %.18e\n",outercount,rk);
   
   if(((rk <= eps) && (rel_prec == 0)) || ((rk <= eps*sourcesquarenorm) && (rel_prec == 1)))
   {
     printf("Reached solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);

     printf("EO Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);
   
  
     stop = clock();
     timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
        
     finalize_mixedsolve();
     return(iter*max_innersolver_it);  // MAYBE ONE SHOULD KEEP TRACK OF REAL INNER SOLVER STEPS
   }
   
  //initialize spin fields on device
  convert2REAL4_spin(g_spinor_field[DUM_SOLVER],h2d_spin);
  
  cudaMemcpy(dev_spinin, h2d_spin, dev_spinsize, cudaMemcpyHostToDevice);
  printf("%s\n", cudaGetErrorString(cudaGetLastError()));

   // solve in single prec on device
   // D p_k = r_k
   printf("Entering inner solver\n");
   assert((startinner = clock())!=-1);
   dev_cg_eo(dev_gf, dev_spinin, dev_spinout, dev_spin1, dev_spin2, dev_spin3, dev_spin4, dev_spin5, dev_grid,dev_nn, dev_output,NULL, T, LZ,0,(REAL) finaleps);
   stopinner = clock();
   timeelapsed = (double) (stopinner-startinner)/CLOCKS_PER_SEC;
   printf("Inner solver done\nTime elapsed: %.6e sec\n", timeelapsed);
 
   // copy back
   cudaMemcpy(h2d_spin, dev_spinout, dev_spinsize, cudaMemcpyDeviceToHost);
   printf("%s\n", cudaGetErrorString(cudaGetLastError()));
   
   convert2double_spin(h2d_spin, g_spinor_field[DUM_SOLVER+2]);
   // x_(k+1) = x_k + p_k
   add(g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+1],g_spinor_field[DUM_SOLVER+2],N);

   outercount ++;   
}// outer loop 
    
     printf("Did NOT reach solver precision of eps=%.2e\n",eps);
     //multiply with Qtm_minus_psi (for non gpu done in invert_eo.c)
     Qtm_minus_psi(g_spinor_field[DUM_SOLVER+3], g_spinor_field[DUM_SOLVER+1]);
     assign(P, g_spinor_field[DUM_SOLVER+3], N);
    

    assert((stop = clock())!=-1);
    timeelapsed = (double) (stop-start)/CLOCKS_PER_SEC;
    printf("Inversion done in mixed precision.\n Number of iterations in outer solver: %d\n Squared residue: %.8e\n Time elapsed: %.6e sec\n", outercount, rk, timeelapsed);

  finalize_mixedsolve();
  return(-1);
}







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






// host/ device interaction

void to_device (dev_spinor * device, spinor * host, dev_spinor * auxiliary, int size) {

  convert2REAL4_spin(host, auxiliary);						// auxiliary = (float) host
  cudaMemcpy(device, auxiliary, size, cudaMemcpyHostToDevice);			// device = auxiliary  (on device)

}


void to_host (spinor * host, dev_spinor * device, dev_spinor * auxiliary, int size) {

  cudaMemcpy(auxiliary, device, size, cudaMemcpyDeviceToHost);			// auxiliary = device  (on device)
  convert2double_spin(auxiliary, host);						// host = (double) auxiliary

}





// puts the additional nd parameters mubar and epsbar on the device
__global__ void he_cg_init_nd_additional (float param_mubar, float param_epsbar) {

  mubar  = param_mubar;
  epsbar = param_epsbar;

}






// derived from Flo's function  dev_mul_one_pm_imu_inv
//	order of the arguments also like Flo's convention: (spinin, spinout)

// applies (1 +- imubar*gamma5)
// uses shared local memory for manipulation	// really ??	where ??
// one thread per lattice site


__global__ void dev_mul_one_pm_imubar_gamma5 (dev_spinor * sin,
                                              dev_spinor * sout,
                                              REAL sign         ) {
   
  dev_spinor slocal[6];									// dev_spinor = float4		// 6*float4 = 24 floats		// auxiliary for each thread
  
  dev_complex pm_imu = dev_initcomplex(0.0, sign * mubar);				// dev_complex = struct { REAL re; REAL im; }	// pm_imu.re = 0.0
  																	// pm_imu.im = sign * mubar
  int pos = threadIdx.x + blockDim.x*blockIdx.x;
  
  if (pos < dev_VOLUME) {
    dev_skalarmult_gamma5_spinor(&(slocal[0]), pm_imu, &(sin[6*pos]) );			// slocal  =  pm_imu * (gamma5) * sin
    dev_add_spinor_assign(&(slocal[0]), &(sin[6*pos]));					// slocal  =  slocal + sin  =  pm_imu * (gamma5) * sin + sin
    dev_realmult_spinor_assign(&(sout[6*pos]), 1.0, &(slocal[0]) );			// sout    =  slocal
  }
}






// initializes and allocates all quantities for the mixed solver
// more precise:
//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
//	allocates memory for all spinor fields
//	puts the nn- and eoidx-fields on device memory

void init_mixedsolve_eo_nd(su3** gf) {	// gf is the full gauge field
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  cudaError_t cudaerr;		// CUDA errors
  int ndev;			// number of devices
  size_t dev_gfsize;		// size of the gauge-field on device memory
  size_t nnsize;		// size of memory for nn-table
  size_t idxsize;		// size of memory for even/odd-positions
  size_t dev_spinsize;		// size of memory for spinors
  int grid[5];			// array for grid specifications
  float * host_output;		// ??
  
  /*
  #ifdef GF_8
    dev_su3_8 * h2d_gf;
  #else
    dev_su3_2v * h2d_gf;
  #endif
  */
  
  /*
  int * nn;
  int * nn_eo;
  int * nn_oe;
  int * eoidx_even;
  int * eoidx_odd;
  */
  
  
  
  
  // get number of devices
  ndev = find_devices();
  if(ndev == 0){
    fprintf(stderr, "Error: no CUDA devices found. Aborting...\n");
    exit(300);
  }
  // try to set active device to device_num given in input file
  if(device_num < ndev){
    printf("Setting active device to: %d\n", device_num);
    cudaSetDevice(device_num);
  }
  else{
    fprintf(stderr, "Error: There is no CUDA device with No. %d. Aborting...\n",device_num);
    exit(301);
  }
  if((cudaerr=cudaGetLastError())!=cudaSuccess){
    printf("Error in init_mixedsolve_eo_nd(): Could not set active device. Aborting...\n");
    exit(302);
  }
  
  
  
  
  #ifdef USETEXTURE
    printf("Using texture references.\n");
  #else
    printf("NOT using texture references.\n");
  #endif
  
  
  
  
  /////////////////
  // GAUGE FIELD //
  /////////////////
  #ifdef GF_8
    /* allocate 8 floats for gf = 2*4*VOLUME float4's*/
    printf("Using GF 8 reconstruction.\n");			// dev_su3_8 = float4
    dev_gfsize = 2*4*VOLUME * sizeof(dev_su3_8);		// allocates for each lattice site and for 4 directions  2*float4 = 8 floats  = 8 real parameters
  #else
    /* allocate 2 rows of gf = 3*4*VOLUME float4's*/
    printf("Using GF 12 reconstruction.\n");			// dev_su3_2v = float4
    dev_gfsize = 3*4*VOLUME * sizeof(dev_su3_2v); 		// allocates for each lattice site and for 4 directions  3*float4 = 12 floats = 2 rows of complex 3-vectors
  #endif
  
  
  if ( (cudaerr = cudaMalloc((void **) &dev_gf, dev_gfsize)) != cudaSuccess ){	// allocates memory for the gauge field dev_gf on device
    printf("Error in init_mixedsolve_eo_nd(): Memory allocation of gauge field failed. Aborting...\n");
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    exit(200);
  }
  else{
    printf("Allocated gauge field on device.\n");
  }  
  
  
  #ifdef GF_8
    // Allocate REAL conversion gf on host ??
    h2d_gf = (dev_su3_8 *)malloc(dev_gfsize); 			// allocates on host
    su3to8(gf,h2d_gf);						// h2d_gf  is the gauge field  gf  with the 8-real-parameter-representation (according to M. Clark, p. 28)
  #else
    // Allocate REAL conversion gf on host ??
    h2d_gf = (dev_su3_2v *)malloc(dev_gfsize);			// allocates on host
    su3to2vf4(gf,h2d_gf);					// h2d_gf  is the gauge field  gf  with the first two rows stored
  #endif
  
  cudaMemcpy(dev_gf, h2d_gf, dev_gfsize, cudaMemcpyHostToDevice);
  								// dev_gf = h2d_gf  on device memory
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Copying dev_gf to device failed.", "Allocated dev_gf on device.");
  		#endif
  
  
  
  
  //////////
  // GRID //
  //////////
  nnsize = 8*VOLUME*sizeof(int);				// size of memory for 8*VOLUME integers
  nn = (int *) malloc(nnsize);					// allocate this memory on host
  nn_eo = (int *) malloc(nnsize/2);				// allocate half this memory
  nn_oe = (int *) malloc(nnsize/2);				// allocate half this memory
  cudaMalloc((void **) &dev_nn, nnsize);			// memory on device
  cudaMalloc((void **) &dev_nn_eo, nnsize/2);			// half the memory on device
  cudaMalloc((void **) &dev_nn_oe, nnsize/2);			// half the memory on device
  
  
  idxsize = VOLUME/2*sizeof(int);				// size of memory necessary for VOLUME/2 integers
  eoidx_even = (int *) malloc(idxsize);				// allocate on host
  eoidx_odd = (int *) malloc(idxsize);				// allocate on host
  cudaMalloc((void **) &dev_eoidx_even, idxsize);		// allocate on device
  cudaMalloc((void **) &dev_eoidx_odd, idxsize);		// allocate on device
  
  
  initnn();							// initialize nearest-neighbour table for gpu
  initnn_eo();							// initialize nearest-neighbour table for gpu with even-odd enabled
  
  
  cudaMemcpy(dev_nn, nn, nnsize, cudaMemcpyHostToDevice);	// copies the previous initialized index-arrays from host to device memory
  cudaMemcpy(dev_nn_eo, nn_eo, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nn_oe, nn_oe, nnsize/2, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_even, eoidx_even, idxsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_eoidx_odd, eoidx_odd, idxsize, cudaMemcpyHostToDevice);
  
  
  
  free(eoidx_odd);						// deallocates the host memory for the field
  free(eoidx_even);						// they are only on the device
  free(nn_oe);
  free(nn_eo);							// not necessary for locally defined variables ??
  free(nn); 
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid stuff failed.", "Allocated grid stuff on device.");
  		#endif
  
  
  
  
  /////////////
  // SPINORS //							// allocates device memory for the odd part of the spinor fields (24 = 4*6 floats per odd lattice sites)
  /////////////							// now we have to consider 2 flavors: up, dn
  
  dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);			// remember: dev_spinor = float4
  
  cudaMalloc((void **) &dev_spin1_up, dev_spinsize);   		// allocates device memory for the fields spinor fields used in dev_cg_eo_nd(...)
  cudaMalloc((void **) &dev_spin1_dn, dev_spinsize);		// pointing to device
  cudaMalloc((void **) &dev_spin2_up, dev_spinsize);		// ...
  cudaMalloc((void **) &dev_spin2_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin3_dn, dev_spinsize);
  /*
  cudaMalloc((void **) &dev_spin4_up, dev_spinsize);		// not needed
  cudaMalloc((void **) &dev_spin4_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin5_dn, dev_spinsize);
  */
  
  cudaMalloc((void **) &dev_spinin_up , dev_spinsize);		// host/device interaction
  cudaMalloc((void **) &dev_spinin_dn , dev_spinsize);		// inner/outer interaction
  cudaMalloc((void **) &dev_spinout_up, dev_spinsize);
  cudaMalloc((void **) &dev_spinout_dn, dev_spinsize);
  
  		// debug	// host code
  		if ( (void *) (h2d_spin_up = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_up. Aborting...\n");
  		  exit(200);
  		}
  		if ( (void *) (h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize) ) == NULL) {
  		  printf("Could not allocate memory for h2d_spin_dn. Aborting...\n");
  		  exit(200);
  		}
  // h2d_spin_up = (dev_spinor *) malloc(dev_spinsize);		// for transferring the spin field in double precision on host to single precision on device
  // h2d_spin_dn = (dev_spinor *) malloc(dev_spinsize);		// on host pointing to host
  
  cudaMalloc((void **) &dev_spin_eo1_up, dev_spinsize);		// used for matrix_multiplication32(...)
  cudaMalloc((void **) &dev_spin_eo1_dn, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_up, dev_spinsize);
  cudaMalloc((void **) &dev_spin_eo3_dn, dev_spinsize);
  /*
  cudaMalloc((void **) &dev_spin_eo2_up, dev_spinsize);		// no memory allocation needed
  cudaMalloc((void **) &dev_spin_eo2_dn, dev_spinsize);		// will point to already allocated memory when used in matrix_multiplication
  */
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of spinor fields failed.", "Allocated spinor fields on device.");
  		#endif
  
  
  
  
  ////////////
  // output //						// ??
  ////////////
  output_size = LZ*T*sizeof(float); 			// parallel in t and z direction
  cudaMalloc((void **) &dev_output, output_size);	// output array
  host_output = (float *) malloc(output_size);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation output stuff failed.", "Allocated output stuff on device.");
  		#endif
  
  
  
  
  ////////////////////////////
  // grid[ ] specifications //							// allocate and initializes the array grid[5] on device
  ////////////////////////////
  grid[0]=LX; grid[1]=LY; grid[2]=LZ; grid[3]=T; grid[4]=VOLUME/2;		// it contains the dimensions of the lattice and the volume of the eo-sublattice
  										// dev_VOLUME is half of VOLUME for eo
 
  cudaMalloc((void **) &dev_grid, 5*sizeof(int));				// dev_grid
  cudaMemcpy(dev_grid, &(grid[0]), 5*sizeof(int), cudaMemcpyHostToDevice);
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in init_mixedsolve_eo_nd(). Memory allocation of grid[] specifications failed.", "Allocated grid[] specifications on device.");
  		#endif
  
  
}//init_mixedsolve_eo_nd()





// deallocates the previous allocated memory

void finalize_mixedsolve_eo_nd(void) {
  
  cudaError_t cudaerr;
  
  cudaFree(dev_spin1_up);
  cudaFree(dev_spin1_dn);
  cudaFree(dev_spin2_up);
  cudaFree(dev_spin2_dn);
  cudaFree(dev_spin3_up);
  cudaFree(dev_spin3_dn);
  /*
  cudaFree(dev_spin4_up);
  cudaFree(dev_spin4_dn);
  cudaFree(dev_spin5_up);
  cudaFree(dev_spin5_dn);
  */
  
  cudaFree(dev_spinin_up);
  cudaFree(dev_spinin_dn);
  cudaFree(dev_spinout_up);
  cudaFree(dev_spinout_dn);
  
  free(h2d_spin_up);
  free(h2d_spin_dn);
  
  cudaFree(dev_spin_eo1_up);
  cudaFree(dev_spin_eo1_dn);
  cudaFree(dev_spin_eo3_up);
  cudaFree(dev_spin_eo3_dn);
  /*
  cudaFree(dev_spin_eo2_up);
  cudaFree(dev_spin_eo2_dn);
  */
  
  
  cudaFree(dev_nn);
  cudaFree(dev_nn_eo);
  cudaFree(dev_nn_oe);
  cudaFree(dev_eoidx_even);
  cudaFree(dev_eoidx_odd);
  

  cudaFree(dev_gf);
  cudaFree(dev_output);
  cudaFree(dev_grid);
  
  
  free(h2d_gf);
  
  
  // Clean up CUDA API for calling thread	// ??
  cudaThreadExit();				// is essential
  
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in finalize_mixedsolve_eo_nd(). Device memory deallocation failed", "Device memory deallocated.");
  		#endif
  
}






///////////////////////////
// MATRIX MULTIPLICATION //
///////////////////////////

// the GPU implementation of  Q_Qdagger_ND(...)  from Nondegenerate_Matrix.c
//	Flo's equivalent function for the standard and non-nd case is  dev_Qtm_pm_psi

void matrix_multiplication32 (dev_spinor * spinout_up, dev_spinor * spinout_dn,
                              dev_spinor * spinin_up , dev_spinor * spinin_dn ,
                              int gridsize1, int blocksize1, int gridsize2, int blocksize2,
                              int gridsize3, int blocksize3, int gridsize4, int blocksize4) {
  
  
  // we will use the auxiliary fields  dev_spin_eo{1,2}_up/dn  for working on and buffering
  // and set  dev_spin_eo2_up/dn  equal  spinout_up/dn
  // spinin_up/dn  have to remain unchanged !!
  // spinout_up/dn  can be freely used
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  int N_sites  =    VOLUME/2;		// #lattice sites
  int N_floats = 24*VOLUME/2;		// #floats
  
  
  
  
  ////////////////////////////////////
  //   MATCHING with Q_Qdagger_ND   //
  ////////////////////////////////////
  //                                //
  // _strange = _up                 //
  // _charm   = _dn                 //
  //                                //
  // DUM_MATRIX   = dev_spin_eo1_up //
  // DUM_MATRIX+1 = dev_spin_eo1_dn //
  //                                //
  // DUM_MATRIX+2 = dev_spin_eo2_up //
  // DUM_MATRIX+3 = dev_spin_eo2_dn //
  //                                //
  ////////////////////////////////////
  
  
  
  
  ///////////////////////////////////
  // INITIALIZATIONS & ASSIGNMENTS //	// have to use (one) other auxiliary field(s) than the calling function dev_cg_eo_nd
  ///////////////////////////////////
  
  dev_spin_eo2_up = spinout_up;		// need no memory allocated
  dev_spin_eo2_dn = spinout_dn;
  												///////////// THEORY ////////////////////////////////////////////////////////////////
  												//                                                                                 //
  												//  (Q_tilde) = gamma5 * ((M_oo) - (M_oe)(Mee^-1)(M_eo))                           //
  												//  (Q_tilde)(Q_tilde_dagger) * (up,dn) = (Q_tilde) * (b,a)                        //
  ///////////////										//                                                    (a,b) = (Q_tilde) * (dn,up)  //
  // MAIN BODY //										//                                                                                 //
  ///////////////										/////////////////////////////////////////////////////////////////////////////////////
  
  
  double nrm = 1.0 / (1.0 + g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  
  ///////////////////////////////////////							/////////////////////////////////
  //        Q_tilde_dagger(2x2)        //							// (a,b) = (Q_tilde) * (dn,up) //
  ///////////////////////////////////////							/////////////////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(spinin_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_dn, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * spinin_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(spinin_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, spinin_up, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * spinin_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);		// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);		// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * spinin_up
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * spinin_dn  +  epsbar * (M_eo) * spinin_up
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * spinin_up  +  epsbar * (M_eo) * spinin_dn
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * spinin_dn  +  nrm*epsbar*(M_eo) * spinin_up
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * spinin_up  +  nrm*epsbar*(M_eo) * spinin_dn
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  #endif									
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  +  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  #endif
  
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_dn, dev_spin_eo2_up, +1.0);			// dev_spin_eo2_up = (1 + imubar) * spinin_dn
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(spinin_up, dev_spin_eo2_dn, -1.0);			// dev_spin_eo2_dn = (1 - imubar) * spinin_up
  
  
  // CUBLAS:													// remember: this is (M_oo) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_up, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*spinin_up  =  (1+imubar)*spinin_dn - epsbar*spinin_up
  cublasSaxpy (N_floats, -g_epsbar, (float *) spinin_dn, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*spinin_dn  =  (1-imubar)*spinin_up - epsbar*spinin_dn
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:													// this is ((M_oo) - (M_oe)(Mee^-1)(M_eo)) * (spinin_dn, spinin_up):
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * spinin_dn  -                    epsbar * spinin_up
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * spinin_dn  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_up
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * spinin_up  -                    epsbar * spinin_dn
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * spinin_up  -  (M_oe)*nrm*epsbar*(M_eo) * spinin_dn
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  	
  
  
  
  
  
  ////////////////////
  // (a,b) -> (b,a) //
  ////////////////////
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_dn, dev_spin_eo3_up);			// dev_spin_eo3_up = dev_spin_eo2_dn
  dev_copy_spinor_field<<<gridsize4, blocksize4>>>(dev_spin_eo2_up, dev_spin_eo3_dn);			// dev_spin_eo3_dn = dev_spin_eo2_up
  
  
  
  
  
  
  ///////////////////////////////////								///////////////////////
  //        Q_tilde(2x2)           //								// (Q_tilde) * (b,a) //
  ///////////////////////////////////								///////////////////////
  
  
  ////////////////////////////
  // (M_oe) (Mee^-1) (M_eo) //
  ////////////////////////////
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_up,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_up, dev_spin_eo1_up, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_up = (M_eo) * dev_spin_eo3_up
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo3_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo3_dn, dev_spin_eo1_dn, dev_eoidx_even, dev_eoidx_odd, dev_nn_eo, 0);	// dev_spin_eo1_dn = (M_eo) * dev_spin_eo3_dn
  #ifdef USETEXTURE
    unbind_texture_spin(1);
  #endif
  
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_up, dev_spin_eo2_up, -1.0);			// dev_spin_eo2_up  =  (1 - imubar) * dev_spin_eo1_up  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo1_dn, dev_spin_eo2_dn, +1.0);			// dev_spin_eo2_dn  =  (1 + imubar) * dev_spin_eo1_dn  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn
  
  
  // CUBLAS:
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_up, 1);
  										// dev_spin_eo2_up  =  dev_spin_eo2_up  +  epsbar * dev_spin_eo1_dn  =  (1 - imubar)*(M_eo) * dev_spin_eo3_up  +  epsbar * (M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, g_epsbar, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_dn, 1);
  										// dev_spin_eo2_dn  =  dev_spin_eo2_dn  +  epsbar * dev_spin_eo1_up  =  (1 + imubar)*(M_eo) * dev_spin_eo3_dn  +  epsbar * (M_eo) * dev_spin_eo3_up
  
  // CUBLAS:
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_up, 1);			// dev_spin_eo2_up  =  nrm * dev_spin_eo2_up  =  nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  
  cublasSscal (N_floats, nrm, (float *) dev_spin_eo2_dn, 1);			// dev_spin_eo2_dn  =  nrm * dev_spin_eo2_dn  =  nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  // Flo:
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_up,1);									// remember: this is ((M_oe) (Mee^-1) (M_eo)) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_up, dev_spin_eo1_up, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_up  =  (M_oe) * dev_spin_eo2_up  =  (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  #endif
  
  #ifdef USETEXTURE
    bind_texture_spin(dev_spin_eo2_dn,1);
  #endif
    dev_Hopping_Matrix<<<gridsize1, blocksize1>>>(dev_gf, dev_spin_eo2_dn, dev_spin_eo1_dn, dev_eoidx_odd, dev_eoidx_even, dev_nn_oe, 1);
  #ifdef USETEXTURE
    unbind_texture_spin(1);							// dev_spin_eo1_dn  =  (M_oe) * dev_spin_eo2_dn  =  (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  +  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  #endif
  
  
  
  ////////////
  // (M_oo) //
  ////////////
  
  // written:
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_up, dev_spin_eo2_up, +1.0);				// dev_spin_eo2_up = (1 + imubar) * dev_spin_eo3_up
  dev_mul_one_pm_imubar_gamma5<<<gridsize2, blocksize2>>>(dev_spin_eo3_dn, dev_spin_eo2_dn, -1.0);				// dev_spin_eo2_dn = (1 - imubar) * dev_spin_eo3_dn
  
  
  // CUBLAS:											// remember: this is (M_oo) * (dev_spin_eo3_up, dev_spin_eo3_dn):
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_dn, 1, (float *) dev_spin_eo2_up, 1);
  												// dev_spin_eo2_up  =  dev_spin_eo2_up - epsbar*dev_spin_eo3_dn  =  (1+imubar)*dev_spin_eo3_up - epsbar*dev_spin_eo3_dn
  cublasSaxpy (N_floats, -g_epsbar, (float *) dev_spin_eo3_up, 1, (float *) dev_spin_eo2_dn, 1);
  												// dev_spin_eo2_dn  =  dev_spin_eo2_dn - epsbar*dev_spin_eo3_up  =  (1-imubar)*dev_spin_eo3_dn - epsbar*dev_spin_eo3_up
  
  
  
  ///////////////////////////////////////
  // (M_oo)  -  (M_oe) (Mee^-1) (M_eo) //
  ///////////////////////////////////////
  
  // CUBLAS:											// this is ( (M_oo) - (M_oe) (Mee^-1) (M_eo) ) * (dev_spin_eo3_up, dev_spin_eo3_dn)
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_up, 1, (float *) dev_spin_eo2_up, 1);
  							// dev_spin_eo2_up  =  dev_spin_eo2_up  -  dev_spin_eo1_up  =                      (1+imubar) * dev_spin_eo3_up  -                    epsbar * dev_spin_eo3_dn
  							//                                                             - (M_oe)*nrm*(1-imubar)*(M_eo) * dev_spin_eo3_up  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_dn
  cublasSaxpy (N_floats, -1.0, (float *) dev_spin_eo1_dn, 1, (float *) dev_spin_eo2_dn, 1);
  							// dev_spin_eo2_dn  =  dev_spin_eo2_dn  -  dev_spin_eo1_dn  =                      (1-imubar) * dev_spin_eo3_dn  -                    epsbar * dev_spin_eo3_up
  							//                                                             - (M_oe)*nrm*(1+imubar)*(M_eo) * dev_spin_eo3_dn  -  (M_oe)*nrm*epsbar*(M_eo) * dev_spin_eo3_up
  
  
  ////////////
  // gamma5 //
  ////////////
  
  // Flo:
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_up, dev_spin_eo2_up);					// dev_spin_eo2_up = gamma5 * dev_spin_eo2_up 
  dev_gamma5<<<gridsize3, blocksize3>>>(dev_spin_eo2_dn, dev_spin_eo2_dn);					// dev_spin_eo2_dn = gamma5 * dev_spin_eo2_dn
  
  
  
  
  /*
  ////////////
  // output //										// output is already done by setting  dev_spin_eo2_up/dn = spinout_up/dn
  ////////////
  
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_up, spinout_up);		// spinout_up = dev_spin_eo2_up
  dev_copy_spinor_field<<<griddim2, blockdim2 >>>(dev_spin_eo2_dn, spinout_dn);		// spinout_dn = dev_spin_eo2_dn
  */
  
  
  return;
  
}//matrix_multiplication32()






// will be used to count the floating point operations per (24 * (#lattice sites)) = (#floats)
// i.e total has to be multiplied by N_floats

void flopcount(unsigned long long int& total, int add) {

  total = total + add;
  
}







extern "C" void benchmark_eo_nd (spinor * const Q_up, spinor * const Q_dn, int N) {

  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //													//
  // total FLOPS  =  (#iterations) * (FLOPS/matrix application) * (#lattice sites)			//
  //													//
  //													//
  // FLOPS per lattice site and application of the function,						//
  // count the floating point op's on device:								//
  //													//
  // dev_Hopping_Matrix	          = 4136								//
  // dev_mul_one_pm_imubar_gamma5 = 120									//
  // dev_gamma5                   = 12									//
  //													//
  // cublasSaxpy                  = 24*2 = 48								//
  // cublasSscal                  = 24*1 = 24								//
  //													//
  //													//
  // (FLOPS/matrix application)  =  2 * (4*4136 + 4*120 + 6*48 + 2*24 + 2*12)  =  2 * 17384  =  34768	//
  //													//
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  double timeelapsed;
  clock_t start_benchmark, stop_benchmark;
  int i;
  double flops;
  double device_flops;
  int flops_per_app = 34768;
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  size_t dev_spinsize = 6*VOLUME/2 * sizeof(dev_spinor);
  
  
  // init_mixedsolve_eo_nd(g_gauge_field);		// only when externally called
  
  
  dev_spinor * A_up;
  dev_spinor * A_dn;
  dev_spinor * B_up;
  dev_spinor * B_dn;
  dev_spinor * C_up;
  dev_spinor * C_dn;
  
  cudaMalloc((void **) &A_up, dev_spinsize);
  cudaMalloc((void **) &A_dn, dev_spinsize);
  cudaMalloc((void **) &B_up, dev_spinsize);
  cudaMalloc((void **) &B_dn, dev_spinsize);
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in benchmark_eo_nd(). Memory allocation of spinor fields failed.");
  		#endif
  
  
  /*
  #ifdef USETEXTURE
    bind_texture_gf(dev_gf);
  #endif
  */
  
  
  /*		// only when externally called
  //Initialize some stuff
  dev_complex h0,h1,h2,h3,mh0, mh1, mh2, mh3;
  
  h0.re = (REAL)ka0.re;    h0.im = -(REAL)ka0.im;
  h1.re = (REAL)ka1.re;    h1.im = -(REAL)ka1.im;
  h2.re = (REAL)ka2.re;    h2.im = -(REAL)ka2.im;
  h3.re = (REAL)ka3.re;    h3.im = -(REAL)ka3.im;
  
  mh0.re = -(REAL)ka0.re;    mh0.im = (REAL)ka0.im;
  mh1.re = -(REAL)ka1.re;    mh1.im = (REAL)ka1.im;
  mh2.re = -(REAL)ka2.re;    mh2.im = (REAL)ka2.im;
  mh3.re = -(REAL)ka3.re;    mh3.im = (REAL)ka3.im;
  
  // try using constant mem for kappas
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex));
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex)); 
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex));
  */
  
  
  
  int blocksize;		// auxiliary
  
  blocksize = 128;
  int blockdim1, griddim1;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim1 = blocksize;
    griddim1  = VOLUME/2/blocksize;
  }
  else {
    blockdim1 = blocksize;
    griddim1  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = 128;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = 128;
  int blockdim3, griddim3;					// passed:	dev_mul_one_pm_imubar_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim3 = blocksize;
    griddim3  = VOLUME/2/blocksize;
  }
  else {
    blockdim3 = blocksize;
    griddim3  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = 128;
  int blockdim4, griddim4;					// passed:	dev_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim4 = blocksize;
    griddim4  = VOLUME/2/blocksize;
  }
  else {
    blockdim4 = blocksize;
    griddim4  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = 128;
  int blockdim5, griddim5;					// passed:	dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim5 = blocksize;
    griddim5  = VOLUME/2/blocksize;
  }
  else {
    blockdim5 = blocksize;
    griddim5  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  
  		//debug
  		printf("Starting a little BENCHMARK. benchmark_eo_nd().\n");
  
  
  
  
  /*		// only when externally called
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK_NO_SUCCESS_MSG("Kernel error in he_cg_init(). Couldn't initialize some stuff.");
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK_NO_SUCCESS_MSG("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.");
  		#endif
  */
  
  
  
  		/*
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK_NO_SUCCESS_MSG(cublasInit(), "CUBLAS error in benchmark_eo_nd(). Couldn't initialize CUBLAS.");
  		#else
  		  cublasInit();
  		#endif
  		*/
  
  
  		// debug
  		printf("Applying the eo-preconditioned matrix %i times.\n", N);
  
  
  to_device(B_up, Q_up, h2d_spin_up, dev_spinsize);
  to_device(B_dn, Q_dn, h2d_spin_dn, dev_spinsize);
  
  
  // timer
  start_benchmark = clock();		// assert((start_benchmark = clock())!=-1);
  
  
  
  for (i = 0; i < N; i++) {
  
    matrix_multiplication32(A_up, A_dn,
                            B_up, B_dn,
                            griddim2, blockdim2,
                            griddim3, blockdim3,
                            griddim4, blockdim4,
                            griddim5, blockdim5);
  
    C_up = B_up;
    C_dn = B_dn;
    B_up = A_up;
    B_dn = A_dn;
    A_up = C_up;
    A_dn = C_dn;
    
  }
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  // CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
  
  
  
  // timer
  stop_benchmark = clock();		// assert((stop_benchmark = clock())!=-1);
  
  
  timeelapsed = double(stop_benchmark - start_benchmark) / double(CLOCKS_PER_SEC);
  device_flops = N * VOLUME/2 * flops_per_app;
  flops = N * VOLUME/2 / timeelapsed / 1.0e9 * flops_per_app;
  
  
  printf("\ttime:        %.2e sec\n", timeelapsed);
  printf("\tflop's:      %.2e flops\n", device_flops);
  printf("\tperformance: %.2e Gflop/s\n", flops);
  
  
  cudaFree(A_up);
  cudaFree(A_dn);
  cudaFree(B_up);
  cudaFree(B_dn);
  
  // finalize_mixedsolve_eo_nd();		// only when externally called
  
  /*
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK_NO_SUCCESS_MSG(cublasShutdown(), "CUBLAS error in benchmark_eo_nd(). Couldn't shut down CUBLAS.");
  		#else
  		  cublasShutdown();
  		#endif
  */
  
  
}//benchmark_eo_nd()










////////////////////////
// CONJUGATE GRADIENT //
////////////////////////

// for the odd field after even/odd-preconditioning
// single precision on GPU

int cg_eo_nd (dev_su3_2v * gf,
              dev_spinor * P_up, dev_spinor * P_dn,
              dev_spinor * Q_up, dev_spinor * Q_dn,
              int max_iter,
              int check_abs , int check_rel,
              double eps_abs, double eps_rel       ) {
  
  // P_up/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  float rr_up;
  float rr_dn;
  float rr;
  float rr_old;
  float r0r0;
  
  float dAd_up;
  float dAd_dn;
  float dAd;
  
  float alpha;
  float beta;
  
  // (auxiliary) device fields
  dev_spinor *  r_up, *  r_dn,
             * Ad_up, * Ad_dn,
             *  x_up, *  x_dn,
             *  d_up, *  d_dn,
             * Ax_up, * Ax_dn;		// for recalculating the residue
  
  // counting
  int j;				// iteration counter
  
  // formal parameters
  size_t dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);
  int N_sites  =    VOLUME/2;
  int N_floats = 24*VOLUME/2;		// (single precision) CUBLAS functions get the number of floats as input
  
  // algorithm control parameters
  int N_recalc_res = 10;		// recalculate residue r(k+1) = b - A*x(k+1) each N_recalc_res iteration
  
  
  
  
  /////////////////////////////////////////////
  // CUDA block- and gridsize specifications //
  /////////////////////////////////////////////
  
  // int gridsize;		// auxiliary
  int blocksize;		// auxiliary
  
  blocksize = BLOCKSIZE1;
  int blockdim1, griddim1;					// here:	dev_zero_spinor_field , dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim1 = blocksize;
    griddim1  = VOLUME/2/blocksize;
  }
  else {
    blockdim1 = blocksize;
    griddim1  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE2;
  int blockdim2, griddim2;					// passed:	dev_Hopping_Matrix
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim2 = blocksize;
    griddim2  = VOLUME/2/blocksize;
  }
  else {
    blockdim2 = blocksize;
    griddim2  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE3;
  int blockdim3, griddim3;					// passed:	dev_mul_one_pm_imubar_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim3 = blocksize;
    griddim3  = VOLUME/2/blocksize;
  }
  else {
    blockdim3 = blocksize;
    griddim3  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE4;
  int blockdim4, griddim4;					// passed:	dev_gamma5
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim4 = blocksize;
    griddim4  = VOLUME/2/blocksize;
  }
  else {
    blockdim4 = blocksize;
    griddim4  = (int) ((VOLUME/2/blocksize) + 1);
  }
  
  blocksize = BLOCKSIZE5;
  int blockdim5, griddim5;					// passed:	dev_copy_spinor_field
  if ( (VOLUME/2) % blocksize == 0 ) {
    blockdim5 = blocksize;
    griddim5  = VOLUME/2/blocksize;
  }
  else {
    blockdim5 = blocksize;
    griddim5  = (int) ((VOLUME/2/blocksize) + 1);
  }
  		/*
  		// debug
  		printf("griddim1 = %i, blockdim1 = %i\n", griddim1, blockdim1);
  		printf("griddim2 = %i, blockdim2 = %i\n", griddim2, blockdim2);
  		printf("griddim3 = %i, blockdim3 = %i\n", griddim3, blockdim3);
  		printf("griddim4 = %i, blockdim4 = %i\n", griddim4, blockdim4);
  		printf("griddim5 = %i, blockdim5 = %i\n", griddim5, blockdim5);
  		*/
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  x_up  = P_up;							// can use the output spinors also as auxiliary fields
  x_dn  = P_dn;							//	saves copying the output spinor field
  /*
  r_up  = Q_up;							// could possibly be done if  Q_up/dn  is not used as feedback
  r_dn  = Q_dn;							//	would save one field and one copying the field
  */
  r_up  = dev_spin1_up;						// use these pointers to the allocated space on device memory (allocated by init_mixedsolve_eo_nd)
  r_dn  = dev_spin1_dn;
  d_up  = dev_spin2_up;
  d_dn  = dev_spin2_dn;
  Ad_up = dev_spin3_up;
  Ad_dn = dev_spin3_dn;
  Ax_up = Ad_up;						// works as long as no initial guess vector x(0) is passed to cg_eo_nd()
  Ax_dn = Ad_dn;
  
  
  
  
  /////////////////////
  // INITIALIZATIONS //
  /////////////////////
  
  /*		// relocated to mixedsolve_eo_nd(), before here were:
  		// Initialize some stuff ...
  		// try using constant mem for kappas ...
  */
  
  /*
  // bind texture gf
  #ifdef USETEXTURE						// needed for subfunctions of dev_Hopping_Matrix(...)
    bind_texture_gf(gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  #endif							//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
    		#endif
  */
  
  /*		// relocated to mixedsolve_eo_nd(), before here were:
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  */
  
  /*
  // cublasInit();			// init CUBLAS
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#else
  		  cublasInit();
  		#endif
  */
  
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  
  // initialize x(0) = 0	// will be added up
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(x_up);
  dev_zero_spinor_field<<<griddim1, blockdim1>>>(x_dn);
  
  
  // r(0) = b - A*x(0) = b
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_up, r_up);
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(Q_dn, r_dn);
  
  
  // d(0) = r(0)
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(r_up, d_up);
  dev_copy_spinor_field<<<griddim1, blockdim1>>>(r_dn, d_dn);
  
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in cg_eo_nd(). Initializing spinor fields on device failed.", "Spinor fields initialized on device.");
  		#endif
  
  
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = cublasSdot(N_floats, (float *) r_up, 1, (float *) r_up, 1);
  rr_dn = cublasSdot(N_floats, (float *) r_dn, 1, (float *) r_dn, 1);
  rr    = rr_up + rr_dn;
  
  		// benchmark
  		#ifdef GPU_BENCHMARK
  		  flopcount(device_flops, 2*2);
  		  // flopcount(device_flops, 2*2*N_floats);
  		#endif
  
  
  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  
  //////////
  // LOOP //
  //////////
  
  		// debug
    		printf("\nEntering inner loop.\n");
  
  		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  // CUBLAS_CORE_CHECK("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.", "Initial inner residue calculated.");
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Calculating initial residue failed.");
		#endif
  
  		// debug
  		printf("Initial inner residue: %.6e\n", r0r0);
  
  
  for (j = 0; j < max_iter; j++) {
    
    
    #ifndef MATRIX_DEBUG
    
    // A*d(k)
    matrix_multiplication32(Ad_up, Ad_dn,										// normally:  matrix_multiplication32()
                             d_up,  d_dn,										// debugging: matrix_debug1(), matrix_multiplication_test()
                            griddim2, blockdim2,
                            griddim3, blockdim3,
                            griddim4, blockdim4,
                            griddim5, blockdim5);
    
  		// debug	// CUDA		// also other stuff ?!
  		#ifdef CUDA_DEBUG
  		  // CUDA_CHECK("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.", "The matrix was applied on GPU.");
  		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in matrix_muliplication32(). Applying the matrix on GPU failed.");
  		#endif
    		
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 1448);
      		  // flopcount(device_flops, 1448*N_floats);
    		#endif
    
    #else
    
    		// debug	// apply the host matrix on trial
    		
    		// host/device interaction
    		to_host(g_chi_up_spinor_field[DUM_SOLVER+3], d_up, h2d_spin_up, dev_spinsize);
    		to_host(g_chi_dn_spinor_field[DUM_SOLVER+3], d_dn, h2d_spin_dn, dev_spinsize);
    		
    		// matrix multiplication
    		printf("This is Q_Qdagger_ND(). ");
    		Q_Qdagger_ND(g_chi_up_spinor_field[DUM_SOLVER+4], g_chi_dn_spinor_field[DUM_SOLVER+4],			// normally:  Q_Qdagger_ND()
    		             g_chi_up_spinor_field[DUM_SOLVER+3], g_chi_dn_spinor_field[DUM_SOLVER+3] );		// debugging: matrix_debug2(), Zwitter1(), Zwitter2(), Zwitter3()
    		
    		// host/device interaction
    		to_device(Ad_up, g_chi_up_spinor_field[DUM_SOLVER+4], h2d_spin_up, dev_spinsize);
    		to_device(Ad_dn, g_chi_dn_spinor_field[DUM_SOLVER+4], h2d_spin_dn, dev_spinsize);
    		
    		
    				// debug	// CUDA
  				#ifdef CUDA_DEBUG
  			  	// CUDA_CHECK("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.", "The matrix was applied on CPU.");
  			  	CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in cg_eo_nd(). Applying the matrix on CPU failed.");
  				#endif
    
    #endif
    
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    dAd_up = cublasSdot(N_floats, (float *) d_up, 1, (float *) Ad_up, 1);
    dAd_dn = cublasSdot(N_floats, (float *) d_dn, 1, (float *) Ad_dn, 1);
    dAd    = dAd_up + dAd_dn;
    
    alpha  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    // x(k+1) = x(k) + alpha*d(k)
    cublasSaxpy(N_floats, alpha, (float *) d_up, 1, (float *) x_up, 1);
    cublasSaxpy(N_floats, alpha, (float *) d_dn, 1, (float *) x_dn, 1);
    
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    // r(k+1)
    if ( (j+1) % N_recalc_res != 0 ) {	// r(k+1) = r(k) - alpha*A*d(k)
      cublasSaxpy(N_floats, -1.0*alpha, (float *) Ad_up, 1, (float *) r_up, 1);
      cublasSaxpy(N_floats, -1.0*alpha, (float *) Ad_dn, 1, (float *) r_dn, 1);
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    }
    
    else {				// recalculate residue r(k+1) = b - A*x(k+1)
    					//	"feedback"
      		// debug
      		printf("Recalculating the inner residue.\n");
      // A*x(k+1)
      matrix_multiplication32(Ax_up, Ax_dn,
                               x_up,  x_dn,
                              griddim2, blockdim2,
                              griddim3, blockdim3,
                              griddim4, blockdim4,
                              griddim5, blockdim5);
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 1448);
      		  // flopcount(device_flops, 1448*N_floats);
      		#endif
      
      // r(k+1) = b - A*x(k+1)
      cublasScopy(N_floats, (float *) Q_up, 1, (float *) r_up, 1);		// r_up = Q_up
      cublasScopy(N_floats, (float *) Q_dn, 1, (float *) r_dn, 1);		// r_dn = Q_dn
      cublasSaxpy(N_floats, -1.0, (float *) Ax_up, 1, (float *) r_up, 1);	// r_up = Q_up - Ax_up
      cublasSaxpy(N_floats, -1.0, (float *) Ax_dn, 1, (float *) r_dn, 1);	// r_dn = Q_dn - Ax_dn
      
      		// benchmark
      		#ifdef GPU_BENCHMARK
      		  flopcount(device_flops, 2*2);
      		  // flopcount(device_flops, 2*2*N_floats);
      		#endif
    
    }
        
    
    // r(k+1)*r(k+1)
    rr_up  = cublasSdot(N_floats, (float *) r_up, 1, (float *) r_up, 1);
    rr_dn  = cublasSdot(N_floats, (float *) r_dn, 1, (float *) r_dn, 1);
    rr     = rr_up + rr_dn;
    
		// debug	// CUBLAS core function
		#ifdef CUDA_DEBUG
		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). CUBLAS function failed.");
		#endif
		
		//benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*2);
    		  // flopcount(device_flops, 2*2*N_floats);
    		#endif
    
    
    		// debug
    		printf("inner iteration j = %i: rr = %.6e\n", j, rr);
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in cg_eo_nd(). Inner residue is NaN.\n");
    		  exit(-1);
    		}
    
    
    // aborting ?? // check wether precision is reached ...
    if ( (check_abs)&&(rr <= eps_abs) || (check_rel)&&(rr <= eps_rel*r0r0) ) {
      
      		// debug
      		printf("Finished inner loop because of reached precision.\n");
      
      if ((check_rel)&&(rr <= eps_rel*r0r0)) {
      		// debug
      		printf("Reached relative inner solver precision of eps_rel = %.2e\n", eps_rel);
      }
      if ((check_abs)&&(rr <= eps_abs)) {
      		// debug
      		printf("Reached absolute inner solver precision of eps_abs = %.2e\n", eps_abs);
      }
      
      		//debug
      		printf("Final inner residue: %.6e\n", rr);
      
      /*
      #ifdef USETEXTURE
        unbind_texture_gf();
      #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      
      // cublasShutdown();			// ends CUBLAS
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasInit(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
      
      return(j+1);
    }
    
    
    // beta = r(k+1)*r(k+1) / r(k)*r(k)
    beta = rr / rr_old;
    
    
    rr_old = rr;  // for next iteration
    
    
    // d(k+1) = r(k+1) + beta*d(k)
    cublasSscal (N_floats, beta, (float *) d_up, 1);
    cublasSaxpy (N_floats, 1.0 , (float *) r_up, 1, (float *) d_up, 1);
    
    cublasSscal (N_floats, beta, (float *) d_dn, 1);
    cublasSaxpy (N_floats, 1.0 , (float *) r_dn, 1, (float *) d_dn, 1);
    
    		// debug	// CUBLAS core function
    		#ifdef CUDA_DEBUG
    		  CUBLAS_CORE_CHECK_NO_SUCCESS_MSG("CUBLAS error in cg_eo_nd(). Error in CUBLAS function.");
    		#endif
    		
    		// benchmark
    		#ifdef GPU_BENCHMARK
    		  flopcount(device_flops, 2*3);
    		  // flopcount(device_flops, 2*3*N_floats);
    		#endif
  
  
  }//LOOP
  
  
  		// debug
  		printf("Finished inner loop beacuse of maximal number of inner iterations.\n");
  		printf("Final inner residue: %.6e\n", rr);
  
  /*
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
  
  		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
  
  
  // cublasShutdown();
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasInit(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
  */
  
  return(j+1);
  
}//cg_eo_nd()






//////////////////
// MIXED SOLVER //
//////////////////

// iterative refinement, defect correction
// that function is to replace the call of  cg_her_nd  in  invert_doublet_eo.c
// solves the odd part of the full eo and nd problem
//	more precisely we have to invert  Qhat(2x2)*Qhat(2x2)^dagger
//	multiplying by  Qhat(2x2)^dagger  is done in  invert_doublet_eo.c

extern "C" int mixedsolve_eo_nd (spinor * P_up, spinor * P_dn,
                                 spinor * Q_up, spinor * Q_dn,
                                 int max_iter, double eps_sq, int rel_prec) {
  
  
  // basically  P_up/dn  and  Q_up/dn  could be used as auxiliary fields
  //	P_up/dn  is the output field (and can be used as initial guess)
  //	Q_up/dn  is not used later in the calling  invert_doublet_eo.c
  //		but will be used as feedback in r(k+1) = b - A*x(k+1)
  
  
  		// debug
  		printf("\n\nmixedsolve_eo_nd():\n");
  		
  		printf("SOLVER PARAMETERS:\n");
  		
  		printf("outer:");
  		printf("\tmaximal iterations: %i\n", max_iter);
  		printf("\trelative check?:    %i\n", bool(rel_prec));
  		printf("\tprecision:          %.8e\n", eps_sq);
  		
  		printf("inner:");
  		printf("\tmaximal iterations: %i\n", max_innersolver_it);
  		printf("\tabsolute check?:    %i\n", bool(innersolver_precision_check_abs));
  		printf("\trelative check?:    %i\n", bool(innersolver_precision_check_rel));
  		printf("\tabsolute precision: %.8e\n", innersolver_precision_abs);
  		printf("\trelative precision: %.8e\n", innersolver_precision_rel);
  
  
  
  
  /////////////////////
  // LOCAL VARIABLES //
  /////////////////////
  
  // CUDA
  cudaError_t cudaerr;
  cublasStatus cublasstatus;
  
  // algorithm
  double rr_up;
  double rr_dn;
  double rr;
  double rr_old;
  double r0r0;
  double bb;
  
  // counting
  int i = 0;					// iteration counter
  int innercount;				// latest inner solver iterations
  int outercount = 0;				// total inner solver iterations
  double flops;
  
  // timing
  clock_t startouter, stopouter;
  clock_t startinner, stopinner;
  // double timeelapsed;
  clock_t innerclocks;
  clock_t totalinnerclocks = 0;
  clock_t totalouterclocks = 0;
  
  // (auxiliary) fields
  spinor *  r_up, *  r_dn,
         * Ad_up, * Ad_dn,
         *  x_up, *  x_dn,
         *  d_up, *  d_dn,
         * Ax_up, * Ax_dn;
  
  // formal parameters
  size_t dev_spinsize = 6*VOLUME/2*sizeof(dev_spinor);		// 24 floats per spinor per even lattice site
  int N_sites  =    VOLUME/2;					// Carsten's functions get the number of lattice points as input
  int N_floats = 24*VOLUME/2;
  
  // algorithm control parameters
  bool rbAx = true;						// choose how to calculate r(k+1)
  bool initial_guess = false;					// choose if initial guess
  
  
  
  
  //////////////////
  // INITIALIZING //
  //////////////////
  
  		//debug
  		printf("init_mixedsolve_eo_nd():\n");
  
  init_mixedsolve_eo_nd(g_gauge_field);				// initializes and allocates all quantities for the mixed solver
  								// more precise:
  								//	puts the gauge field on device as "2 rows" or "8 floats" per SU(3)-matrix
  								//	allocates memory for all spinor fields
  								//	puts the nn- and eoidx-fields on device memory
  
  		//debug
  		printf("mixedsolve_eo_nd():\n");
  
  
  // the following initializations are moved from cg_eo_nd():
  
  // Initialize some stuff
  dev_complex h0, h1, h2, h3, mh0, mh1, mh2, mh3;
  
  h0.re  =  (REAL) ka0.re;	h0.im  = -(REAL) ka0.im;	// ka{0-4} are defined in boundary.c
  h1.re  =  (REAL) ka1.re;	h1.im  = -(REAL) ka1.im;	// what is the meaning?
  h2.re  =  (REAL) ka2.re;	h2.im  = -(REAL) ka2.im;
  h3.re  =  (REAL) ka3.re;	h3.im  = -(REAL) ka3.im;
  
  mh0.re = -(REAL) ka0.re;	mh0.im =  (REAL) ka0.im;
  mh1.re = -(REAL) ka1.re;	mh1.im =  (REAL) ka1.im;
  mh2.re = -(REAL) ka2.re;	mh2.im =  (REAL) ka2.im;
  mh3.re = -(REAL) ka3.re;	mh3.im =  (REAL) ka3.im;
  /*
  // try using constant mem for kappas		// constant memory is cached!
  cudaMemcpyToSymbol("dev_k0c", &h0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k1c", &h1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k2c", &h2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_k3c", &h3, sizeof(dev_complex));
  
  cudaMemcpyToSymbol("dev_mk0c", &mh0, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk1c", &mh1, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk2c", &mh2, sizeof(dev_complex));
  cudaMemcpyToSymbol("dev_mk3c", &mh3, sizeof(dev_complex));
  */
  
  
  // bind texture gf
  #ifdef USETEXTURE						// needed for subfunctions of dev_Hopping_Matrix(...)
    bind_texture_gf(dev_gf);					//	e.g. dev_reconstructgf_2vtexref(...), dev_reconstructgf_8texref(...), 
  #endif							//	in general in functions  dev_reconstructgf[...]  with  "tex1Dfetch(gf_tex[...]"
  
  		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  CUDA_CHECK("CUDA error in bind_texture_gf(). Binding GF to texture failed.", "GF bound to texture.");
    		#endif
  
  
  he_cg_init<<< 1, 1 >>> (dev_grid, (REAL) g_kappa, (REAL)(g_mu/(2.0*g_kappa)), h0, h1, h2, h3);
  		// "he" = "host entry"
  		// BEWARE in dev_tm_dirac_kappa we need the true mu (not 2 kappa mu!)	// ??
  		
  		// dev_LX, dev_LY, dev_LZ, dev_T, dev_VOLUME  =  grid[5]  =  dev_grid[5]
  		//	dev_VOLUME  is necessary for many kernel functions as for instance  dev_gamma5()
  		// initializes  mu, kappa and twokappamu  on the device
  		// initializes the strange  dev_k{0-3}, dev_mk{0-3}  as derived from the  ka{0-3}  from boundary.c
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init(). Couldn't initialize some stuff.", "he_cg_init() succeeded.");
  		#endif
  
  		// debug	// check stuff on device
  		#ifdef STUFF_DEBUG
  		  int host_check_LX, host_check_LY, host_check_LZ, host_check_T, host_check_VOLUME;
  		  cudaMemcpyFromSymbol(&host_check_LX, dev_LX, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LY, dev_LY, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_LZ, dev_LZ, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_T, dev_T, sizeof(int));
  		  cudaMemcpyFromSymbol(&host_check_VOLUME, dev_VOLUME, sizeof(int));
  		  // printf("\teven_odd_flag = %i\n", even_odd_flag);
  		  printf("\tOn device:\n");
  		  printf("\tdev_LX = %i\n", host_check_LX);
  		  printf("\tdev_LY = %i\n", host_check_LY);
  		  printf("\tdev_LZ = %i\n", host_check_LZ);
  		  printf("\tdev_T = %i\n", host_check_T);
  		  printf("\tdev_VOLUME = %i/2 ?!= %i\n", host_check_LX*host_check_LY*host_check_LZ*host_check_T, host_check_VOLUME);
  		  
  		  float host_check_mu, host_check_kappa, host_check_twokappamu;
  		  cudaMemcpyFromSymbol(&host_check_mu, mu, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_kappa, kappa, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_twokappamu, twokappamu, sizeof(float));
  		  // printf("\tOn device:\n");
  		  // printf("\tmu = %f\n", host_check_mu);		// not needed for the nd case
  		  printf("\tkappa = %f\n", host_check_kappa);
  		  // printf("\ttwokappamu = %f\n", host_check_twokappamu);
  		#endif
  
  
  he_cg_init_nd_additional<<<1,1>>> (g_mubar, g_epsbar);
  
  		// debug	// kernel
  		#ifdef CUDA_DEBUG
  		  CUDA_KERNEL_CHECK("Kernel error in he_cg_init_nd_additional(). Couldn't initialize some stuff.", "he_cg_init_nd_additional() succeeded.");
  		#endif
  
  		// debug	// check mubar and epsbar on host and device
  		#ifdef STUFF_DEBUG
  		  // printf("\tOn host:\n");
  		  // printf("\tg_mubar = %f\n", g_mubar);
  		  // printf("\tg_epsbar = %f\n", g_epsbar);
  		  
  		  float host_check_mubar, host_check_epsbar;
  		  cudaMemcpyFromSymbol(&host_check_mubar, mubar, sizeof(float));
  		  cudaMemcpyFromSymbol(&host_check_epsbar, epsbar, sizeof(float));
  		  printf("\tOn device:\n");
  		  printf("\tmubar = %f\n", host_check_mubar);
  		  printf("\tepsbar = %f\n", host_check_epsbar);
  		#endif
  
  
  /*		// necessary ??
  // cublasInit();
  
  		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasInit(), "CUBLAS error in cublasInit(). Couldn't initialize CUBLAS.", "CUBLAS initialized.");
  		#else
  		  cublasInit();
  		#endif
  */
  
  
  
  #ifdef OPERATOR_BENCHMARK
    benchmark_eo_nd(Q_up, Q_dn, 1000);
  #endif
  
  
  
  
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  
  x_up = P_up;							// can use the output spinors also as auxiliary fields
  x_dn = P_dn;							//	can use as initial guess at the same time
  
  
  #ifndef CG_DEBUG
  
  r_up  = g_chi_up_spinor_field[DUM_SOLVER];			// use the pre-allocated memory on host memory
  r_dn  = g_chi_dn_spinor_field[DUM_SOLVER];			// allocated by  init_chi_spinor_field.c  and  invert_doublet.c  !?
  d_up  = g_chi_up_spinor_field[DUM_SOLVER+1];			// the fields  g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, ... , +5}]  are used in  cg_her_nd()
  d_dn  = g_chi_dn_spinor_field[DUM_SOLVER+1];
  Ad_up = g_chi_up_spinor_field[DUM_SOLVER+2];
  Ad_dn = g_chi_dn_spinor_field[DUM_SOLVER+2];
  Ax_up = Ad_up;
  Ax_dn = Ad_dn;
  		// debug
  		printf("Now using the fields g_chi_up/dn_spinor_field[DUM_SOLVER{ , +1, +2}] in the mixedsolve_eo_nd().\n");
  
  #else
  
  		r_up  = (spinor *) malloc(24*N_sites*sizeof(double));		// if using cg_her_nd() as the CG, we cannot use the g_chi_up/dn-fields at the same time
  		r_dn  = (spinor *) malloc(24*N_sites*sizeof(double));
  		d_up  = (spinor *) malloc(24*N_sites*sizeof(double));
  		d_dn  = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ad_up = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ad_dn = (spinor *) malloc(24*N_sites*sizeof(double));
  		Ax_up = Ad_up;
  		Ax_dn = Ad_dn;
  				// debug
  				printf("Now allocating new host space for the fields in mixedsolve_eo_nd().\n");
  
  #endif
  
  
  		// benchmark
  		#ifdef GPU_BENCHMARK
  		  device_flops = 0;
  		#endif
  		
  		#ifdef CPU_BENCHMARK
  		  host_flops = 0;
  		#endif
  
  
  
  
  ///////////////
  // ALGORITHM //
  ///////////////
  
  // timer
  startouter = clock();
  
  
  // r(0)
  if (!initial_guess) {		// r(0) = b = Q	// for x(0) = 0
    assign(r_up, Q_up, N_sites);
    assign(r_dn, Q_dn, N_sites);
    printf("x(0) = 0\n");
  }
  else {			// r(0) = b - A*x(0) = Q - A*P
    bb = square_norm(P_up, N_sites, 0) + square_norm(P_dn, N_sites, 0);
    		// benchmark
    		#ifdef CPU_BENCHMARK
    		  flopcount(host_flops, 2*2);
    		  // flopcount(host_flops, 2*2*N_floats);
    		#endif
    printf("bb = %.10e\n", bb);
    if (bb == 0) {
      assign(r_up, Q_up, N_sites);
      assign(r_dn, Q_dn, N_sites);
      printf("x(0) = 0\n");
    }
    else {
      Q_Qdagger_ND(Ax_up, Ax_dn, P_up, P_dn);
      diff(r_up, Q_up, Ax_up, N_sites);
      diff(r_dn, Q_dn, Ax_dn, N_sites);
      		// benchmark
      		#ifdef CPU_BENCHMARK
      		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
      		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
      		#endif
      printf("x(0) != 0\n");
    }
  }
  
  
  // rr = (r_up)^2 + (r_dn)^2
  rr_up = square_norm(r_up, N_sites, 0);
  rr_dn = square_norm(r_dn, N_sites, 0);
  rr = rr_up + rr_dn;
  
  		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2);
  		  // flopcount(host_flops, 2*2*N_floats);
  		#endif
  
  
  r0r0   = rr; // for relative precision
  rr_old = rr; // for the first iteration
  
  		// debug
  		printf("Initial outer residue: %.10e\n", rr_old);
  
  
  // set to zero	// x_up, x_dn  will be added up		// as  x_up/dn = P_up/dn  up to here  P_up/dn  was not changed
  zero_spinor_field(x_up, N_sites);
  zero_spinor_field(x_dn, N_sites);
  
  
  
  
  ////////////////
  // OUTER LOOP //
  ////////////////
  
  		// debug
    		printf("\nEntering outer loop.");
  
  
  do {		// for (i = 0; i < max_iter; i++) {
    
    i++;
  
    		// debug
    		printf("\nouter iteration i = %i\n", i);
    
    
    
    
    #ifndef CG_DEBUG
    
    // host/device interaction
    to_device(dev_spinin_up, r_up, h2d_spin_up, dev_spinsize);
    to_device(dev_spinin_dn, r_dn, h2d_spin_dn, dev_spinsize);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.", "Fields copied to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Host to device interaction failed.");
    		#endif
    
    
    
    
    ////////////////////////////////////
    // INNER LOOP, CONJUGATE GRADIENT //
    ////////////////////////////////////
    
    // timer
    startinner = clock();
    
    		// debug
    		printf("cg_eo_nd():\n");
    
    // solves A*p(k+1) = r(k)
    //        A*p(0)   = r(0) = b
    innercount = cg_eo_nd(dev_gf,
                          dev_spinout_up, dev_spinout_dn,
                          dev_spinin_up , dev_spinin_dn,
                          max_innersolver_it,
                          innersolver_precision_check_abs, innersolver_precision_check_rel,
                          innersolver_precision_abs      , innersolver_precision_rel      );
    
    outercount = outercount + innercount;
    
    // timer
    stopinner = clock();
    innerclocks = stopinner-startinner;
    totalinnerclocks = totalinnerclocks + innerclocks;
    
    		// debug
    		printf("Inner solver done in: %.4e sec\n", double(innerclocks) / double(CLOCKS_PER_SEC));
    
    
    // host/device interaction
    to_host(d_up, dev_spinout_up, h2d_spin_up, dev_spinsize);
    to_host(d_dn, dev_spinout_dn, h2d_spin_dn, dev_spinsize);
    
    		// debug	// CUDA
    		#ifdef CUDA_DEBUG
    		  // CUDA_CHECK("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.", "Fields copied back to device.");
    		  CUDA_CHECK_NO_SUCCESS_MSG("CUDA error in mixedsolve_eo_nd(). Device to host interaction failed.");
    		#endif
    
    
    #else
    
    
    				// debug
    				printf("cg_her_nd():\n");
    		
    		innercount = cg_her_nd(d_up, d_dn, r_up, r_dn,		// MISTAKE, was: r_up, r_dn, d_up, d_dn,
				           1000, eps_sq/2, 0,
				           VOLUME/2, &Q_Qdagger_ND, 0, 1000);
    		
    		outercount = outercount + innercount;
    		
    				// debug
    				printf("cg_her_nd() on host was used for debugging purposes.\n");
    
    
    #endif
    
    
    		// debug
    		printf("mixedsolve_eo_nd():\n");
    
    
    // x(k+1) = x(k) + d(k+1)
    add(x_up, x_up, d_up, N_sites);
    add(x_dn, x_dn, d_dn, N_sites);
    
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2);
  		  // flopcount(host_flops, 2*N_floats);
  		#endif
    
    
    // r(k+1)
    if (rbAx) {				// r(k+1) = b - A*x(k+1)
      // A*x(k+1)
      Q_Qdagger_ND(Ax_up, Ax_dn, x_up, x_dn);
      		// debug
      		printf("The matrix was applied on CPU in double precision. r = b - Ax\n");
      diff(r_up, Q_up, Ax_up, N_sites);
      diff(r_dn, Q_dn, Ax_dn, N_sites);
    }
    else {				// r(k+1) = r(k) - A*d(k+1)	// makes actually no sense ;)
      // A*d(k+1)
      Q_Qdagger_ND(Ad_up, Ad_dn, d_up, d_dn);
    		// debug
    		printf("The matrix was applied on CPU in double precision. r = r - Ad\n");
      // r(k+1) = r(k) - A*d(k+1)
      diff(r_up, r_up, Ad_up, N_sites);
      diff(r_dn, r_dn, Ad_dn, N_sites);
    }
    
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2*(55+2+2+1+55) + 2);
  		  // flopcount(host_flops, 2*2*(55+2+2+1+55)*N_floats + 2*N_floats);
  		#endif
    
    
    // rr = (rr_up)^2 + (r_dn)^2
    rr_up = square_norm(r_up, N_sites, 0);
    rr_dn = square_norm(r_dn, N_sites, 0);
    rr    = rr_up + rr_dn;
    
    		// debug
    		printf("Outer residue in the outer iteration i = %i after %i total inner iterations : %.10e\n", i, outercount, rr);
    		
    		// debug	// is NaN ?
    		if isnan(rr) {
    		  printf("Error in mixedsolve_eo_nd(). Outer residue is NaN.\n");
    		  exit(-1);
    		}
    		
    		// benchmark
  		#ifdef CPU_BENCHMARK
  		  flopcount(host_flops, 2*2);
  		  // flopcount(host_flops, 2*2*N_floats);
  		#endif
    
    
    // aborting ?? // check wether precision is reached ...
    if ( ((rr <= eps_sq) && (rel_prec == 0))  ||  ((rr <= eps_sq*r0r0) && (rel_prec == 1)) ) {
      
      // timer
      stopouter = clock();
      totalouterclocks = stopouter-startouter - totalinnerclocks;
      
      		/*
      		// benchmark
  		#ifdef GPU_BENCHMARK2
  		  device_flops = 0;
  		  int help = ( 4 + outercount*(1448+5*4+6) + outercount/10*1448 ) * N_floats;
  		  flopcount(device_flops, help);			// N_recalcres = 10
  		#endif
  		*/
      
      		// debug
      		printf("\nEO inversion done in mixed precision.\n");
      		if (rel_prec == 0) printf("Finished outer loop because of reached absolute outer solver precision.\n");
      		if (rel_prec == 1) printf("Finished outer loop because of reached relative outer solver precision.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter) / double(CLOCKS_PER_SEC));
      		// benchmark
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats));
      		  printf("\touter solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      
      
      #ifdef USETEXTURE
        unbind_texture_gf();
      #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      /*
      // cublasShutdown();
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasShutdown(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
      
      		// debug
      		printf("finalize_mixedsolve_eo_nd():\n");
      
      finalize_mixedsolve_eo_nd();
      
      		// debug
      		printf("\n");
      
      return(outercount);
    }
    
    
  }//OUTER LOOP
  while (outercount <= max_iter);
  
  
  // multiplying with Qhat(2x2)^dagger is done in invert_doublet_eo.c
  
  
  // timer
  stopouter = clock();
  totalouterclocks = stopouter-startouter - totalinnerclocks;
  
  		/*
  		// benchmark
  		#ifdef GPU_BENCHMARK2
  		  device_flops = 0;
  		  int help = ( 4 + outercount*(1448+5*4+6) + outercount/10*1448 ) * N_floats;
  		  flopcount(device_flops, help);			// N_recalcres = 10
  		#endif
  		*/
  
  		// debug
  		printf("\nEO inversion done in mixed precision.\n");
  		printf("Finished outer loop, because of maximal number of outer iterations.\n");
      		printf("Total number of inner iterations: %i\n", outercount);
      		printf("Total number of outer iterations: %i\n", i+1);
      		printf("Squared residue: %.10e\n", rr); 
      		printf("Outer solver done in: %.4e sec\n", double(stopouter-startouter)/CLOCKS_PER_SEC);
      		// benchmark
      		#if defined(GPU_BENCHMARK) || defined(GPU_BENCHMARK2)
      		  // REMARK: device_flops has to be multiplied by N_floats !!
      		  flops = device_flops * N_floats / (double(totalinnerclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Inner solver BENCHMARK:\n");
      		  printf("\ttotal inner solver time:   %.2e sec\n", double(totalinnerclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(device_flops) * double(N_floats));
      		  printf("\tinner solver performance:  %.2e Gflop/s\n", flops);
      		#endif
      		
      		#ifdef CPU_BENCHMARK
      		  // REMARK: host_flops has to be multiplied by N_floats !!
      		  flops = host_flops * N_floats / (double(totalouterclocks)/double(CLOCKS_PER_SEC)) / 1.0e9;
      		  printf("Outer solver BENCHMARK:\n");
      		  printf("\ttotal outer solver time:   %.2e sec\n", double(totalouterclocks) / double(CLOCKS_PER_SEC));
      		  printf("\tfloating point operations: %.2e flops\n", double(host_flops) * double(N_floats));
      		  printf("\touter solver performance:  %.2e Gflop/s\n", flops);
      		#endif
  
  
  #ifdef USETEXTURE
    unbind_texture_gf();
  #endif
      
      		// debug	// CUDA
  		#ifdef CUDA_DEBUG
  		  CUDA_CHECK("CUDA error in unbind_texture(). Unbindung the GF texture failed.", "GF texture unbound.");
  		#endif
      
      /*
      // cublasShutdown();
      
      		// debug	// CUBLAS helper function
  		#ifdef CUDA_DEBUG
  		  CUBLAS_HELPER_CHECK(cublasShutdown(), "CUBLAS error in cublasShutdown(). Couldn't shut down CUBLAS.", "CUBLAS is shutted down.");
  		#else
  		  cublasShutdown();
  		#endif
      */
  
  		// debug
  		printf("finalize_mixedsolve_eo_nd():\n");  
  
  finalize_mixedsolve_eo_nd();
  
  		// debug
  		printf("\n");
  
  return(outercount);
  
  
}//mixedsolve_eo_nd()










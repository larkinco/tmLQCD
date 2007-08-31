
/*******************************************************************************
 * $Id$
 *
 * Action of a Dirac operator D (Wilson or twisted) on a given spinor field
 *
 * The externally accessible function is
 *
 *   void D_psi(spinor * const P, spinor * const Q, const double m0)
 *     Computes D*Q and stores the result in P
 *
 * original code by
 *       Martin Luescher <luscher@mail.desy.de>
 *
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "sse.h"
#include "boundary.h"
#ifdef MPI

#endif
#include "D_psi.h"

#if (defined SSE2 || defined SSE3)

static spinor rs __attribute__ ((aligned (16)));

/* Serially Checked ! */
void D_psi(spinor * const P, spinor * const Q){
  int ix,iy,iz;
  su3 *up,*um;
  spinor *s,*sp,*sm,*rn;
  static complex fact1, fact2;

  if(P==Q){
    printf("Error in D_psi (operator.c):\n");
    printf("Arguments must be differen spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }


# if defined MPI

# endif

  fact1.re = 1.;
  fact1.im = g_mu;
  fact2.re = 1.;
  fact2.im = -g_mu;

  iy=g_iup[0][0];
  sp=(spinor *) Q + iy;
  up=&g_gauge_field[0][0];
   
  /************************ loop over all lattice sites *************************/
   
  for (ix=0;ix<VOLUME;ix++){
    s=(spinor *) Q + ix;
    _prefetch_spinor(s);

    /******************************* direction +0 *********************************/

    iy=g_idn[ix][0];
      
    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm);       

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_0);
    _sse_store_up(rs.s2);

    _sse_load_up((*s).s0);
    _sse_vector_cmplx_mul(fact1);
/*     _sse_vector_mul(fact1); */
    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load_up((*s).s2);
    _sse_vector_cmplx_mul(fact2);
/*     _sse_vector_mul(fact1);       */
    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);      
      
    um=&g_gauge_field[iy][0];
    _prefetch_su3(um);
      
    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_add();
      
    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_0);
    _sse_store_up(rs.s3);
    
    _sse_load_up((*s).s1);
    _sse_vector_cmplx_mul(fact1);
/*     _sse_vector_mul(fact1); */
    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load_up((*s).s3);
    _sse_vector_cmplx_mul(fact2);
/*     _sse_vector_mul(fact1);       */
    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3); 

    /******************************* direction -0 *********************************/

    iy=g_iup[ix][1];

    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_0);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);
      
    up+=1;
    _prefetch_su3(up);
      
    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_0);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);
      
    /******************************* direction +1 *********************************/

    iy=g_idn[ix][1];
      
    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_1);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s3); 
      
    um=&g_gauge_field[iy][1];
    _prefetch_su3(um);

    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);       

    /******************************* direction -1 *********************************/

    iy=g_iup[ix][2];

    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_1);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);

    up+=1;
    _prefetch_su3(up);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_1);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s2);

    /******************************* direction +2 *********************************/

    iy=g_idn[ix][2];

    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s3);
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_2);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_add();
    _sse_store(rs.s3);
      
    um=&g_gauge_field[iy][2];
    _prefetch_su3(um);

    _sse_load((*sp).s1);
    _sse_load_up((*sp).s2);
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_2);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_sub();
    _sse_store(rs.s2);      

    /******************************* direction -2 *********************************/

    iy=g_iup[ix][3];

    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s3);
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_2);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s3);
    _sse_vector_sub();
    _sse_store(rs.s3);
      
    up+=1;
    _prefetch_su3(up);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s2);
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_2);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s2);
    _sse_vector_add();
    _sse_store(rs.s2);      
      
    /******************************* direction +3 *********************************/

    iy=g_idn[ix][3];

    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm);

    _sse_load((*sp).s0);
    _sse_load_up((*sp).s2);
    _sse_vector_i_mul();
    _sse_vector_add();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_3);
    _sse_load(rs.s0);
    _sse_vector_add();
    _sse_store(rs.s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_sub();
    _sse_store(rs.s2);
      
    um=&g_gauge_field[iy][3];
    _prefetch_su3(um);

    _sse_load((*sp).s1);
    _sse_load_up((*sp).s3);
    _sse_vector_i_mul();
    _sse_vector_sub();

    _sse_su3_multiply((*up));
    _sse_vector_cmplx_mul(phase_3);
    _sse_load(rs.s1);
    _sse_vector_add();
    _sse_store(rs.s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_add();
    _sse_store(rs.s3);
      
    /******************************* direction -3 *********************************/

    iz=(ix+1+VOLUME)%VOLUME;

    iy=g_iup[iz][0];
      
    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp);

    _sse_load((*sm).s0);
    _sse_load_up((*sm).s2);
    _sse_vector_i_mul();
    _sse_vector_sub();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_3);
    rn = (spinor *) P + ix;
      
    _sse_load(rs.s0);
    _sse_vector_add();
/*     _sse_vector_mul(fact2); */
    _sse_store_nt((*rn).s0);

    _sse_load(rs.s2);
    _sse_vector_i_mul();      
    _sse_vector_add();
/*     _sse_vector_mul(fact2);       */
    _sse_store_nt((*rn).s2);

    up=&g_gauge_field[iz][0];
    _prefetch_su3(up);

    _sse_load((*sm).s1);
    _sse_load_up((*sm).s3);
    _sse_vector_i_mul();
    _sse_vector_add();
      
    _sse_su3_inverse_multiply((*um));
    _sse_vector_cmplxcg_mul(phase_3);
    _sse_load(rs.s1);
    _sse_vector_add();
/*     _sse_vector_mul(fact2);       */
    _sse_store_nt((*rn).s1);

    _sse_load(rs.s3);
    _sse_vector_i_mul();      
    _sse_vector_sub();
/*     _sse_vector_mul(fact2); */
    _sse_store_nt((*rn).s3);
      
    /******************************** end of loop *********************************/

  }
}

#elif ((defined BGL) && (defined XLC))

/**********************************
 *
 * Blue Gene/L Version
 *
 * Author: Carsten Urbach
 *
 **********************************/

void D_psi(spinor * const P, spinor * const Q){
  int ix,iy,iz;
  static complex fact1;
  su3 * restrict up ALIGN;
  su3 * restrict um ALIGN;
  spinor * restrict s ALIGN;
  spinor * restrict sp ALIGN;
  spinor * restrict sm ALIGN;
  spinor * restrict rn ALIGN;
  /* We have 32 registers available */
  double _Complex reg00, reg01, reg02, reg03, reg04, reg05;
  double _Complex reg10, reg11, reg12, reg13, reg14, reg15;
  /* For the gauge field, reuse the first three!*/
  double _Complex u00, u01, u02, u10, u11, u12;
  double _Complex reg20, reg21;
  /* The following contains the result spinor (12 regs) */
  double _Complex rs00, rs01, rs02, rs10, rs11, rs12, rs20, rs21, rs22, 
    rs30, rs31, rs32;

#pragma disjoint(*s, *sp, *sm, *rn, *up, *um, *P, *Q)

  __alignx(16,P);
  __alignx(16,Q);

#    if (defined MPI && !(defined _NO_COMM))

#    endif

  fact1.re = 1.;
  fact1.im = g_mu;

  iy=g_iup[0][0];
  sp=(spinor *) Q + iy;
  up=&g_gauge_field[0][0];

  /**************** loop over all lattice sites ******************/
  for(ix = 0; ix < VOLUME; ix++){
    s=(spinor *) Q + ix;
    rn = (spinor *) P + ix;
    /*********************** direction +0 ************************/

    iy=g_idn[ix][0]; 

    um=&g_gauge_field[iy][0]; 

    _prefetch_su3(um); 
    sm = (spinor*) Q + iy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_add_reg1();
    /* result is now in regx0, regx1, regx2 x = 0,1 */

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(phase_0);
    _bgl_load_rs0((*s).s0);
    _bgl_load_rs1((*s).s1);
    _bgl_load_rs2((*s).s2);
    _bgl_load_rs3((*s).s3);
    _bgl_vector_cmplx_mul_rs(fact1);
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs3_reg1();

    /*********************** direction -0 ************************/

    iy=g_iup[ix][1]; 

    up+=1;
    _prefetch_su3(up); 
    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka0);

    _bgl_add_to_rs0_reg0();
    _bgl_sub_from_rs2_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs3_reg1();

    /*********************** direction +1 ************************/

    iy=g_idn[ix][1]; 

    um=&g_gauge_field[iy][1]; 

    _prefetch_su3(um); 
    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_add_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka1);

    _bgl_add_to_rs0_reg0();
    _bgl_i_mul_sub_from_rs3_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg1();

    /*********************** direction -1 ************************/

    iy=g_iup[ix][2]; 

    up+=1;
    _prefetch_su3(up); 
    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_sub_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka1);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_add_to_rs3_reg0();
    _bgl_i_mul_add_to_rs2_reg1();      

    /*********************** direction +2 ************************/

    iy=g_idn[ix][2]; 

    um=&g_gauge_field[iy][2]; 
    _prefetch_su3(um); 
    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg1_up((*sp).s2);
    _bgl_load_reg0_up((*sp).s3);
    _bgl_vector_add_reg0();
    _bgl_vector_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka2);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_sub_from_rs2_reg1();
    _bgl_add_to_rs3_reg0();


    /*********************** direction -2 ************************/

    iy=g_iup[ix][3]; 

    up+=1;
    _prefetch_su3(up); 
    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg1_up((*sm).s2);
    _bgl_load_reg0_up((*sm).s3);
    _bgl_vector_sub_reg0();
    _bgl_vector_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka2);
      
    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_add_to_rs2_reg1();
    _bgl_sub_from_rs3_reg0();

    /*********************** direction +3 ************************/

    iy=g_idn[ix][3]; 

    um=&g_gauge_field[iy][3]; 
    _prefetch_su3(um); 
    sm = (spinor *) Q + iy;
    _prefetch_spinor(sm); 

    _bgl_load_reg0((*sp).s0);
    _bgl_load_reg1((*sp).s1);
    _bgl_load_reg0_up((*sp).s2);
    _bgl_load_reg1_up((*sp).s3);
    _bgl_vector_i_mul_add_reg0();
    _bgl_vector_i_mul_sub_reg1();

    _bgl_su3_multiply_double((*up));
    _bgl_vector_cmplx_mul_double(ka3);

    _bgl_add_to_rs0_reg0();
    _bgl_add_to_rs1_reg1();
    _bgl_i_mul_sub_from_rs2_reg0();
    _bgl_i_mul_add_to_rs3_reg1();

    /*********************** direction -3 ************************/

    iz=(ix+1+VOLUME)%VOLUME;

    iy=g_iup[iz][0];
      
    up=&g_gauge_field[iz][0];
    _prefetch_su3(up); 
    sp = (spinor *) Q + iy;
    _prefetch_spinor(sp); 

    _bgl_load_reg0((*sm).s0);
    _bgl_load_reg1((*sm).s1);
    _bgl_load_reg0_up((*sm).s2);
    _bgl_load_reg1_up((*sm).s3);
    _bgl_vector_i_mul_sub_reg0();
    _bgl_vector_i_mul_add_reg1();
      
    _bgl_su3_inverse_multiply_double((*um));
    _bgl_vector_cmplxcg_mul_double(ka3);
      
    _bgl_add_to_rs0_reg0();
    _bgl_store_rs0((*rn).s0);
    _bgl_i_mul_add_to_rs2_reg0();
    _bgl_store_rs2((*rn).s2);

    _bgl_add_to_rs1_reg1();
    _bgl_store_rs1((*rn).s1);
    _bgl_i_mul_sub_from_rs3_reg1();
    _bgl_store_rs3((*rn).s3);

    /************************ end of loop ************************/
  }
}


#else

/* Serially Checked ! */

static complex rho1,rho2;
static su3_vector psi,chi;

void D_psi(spinor * const P, spinor * const Q){
  int ix,iy;
  static spinor r;
  su3 * restrict up,* restrict um;
  spinor * restrict rr,* restrict s,* restrict sp,* restrict sm;

  if(P==Q){
    printf("Error in D_psi (operator.c):\n");
    printf("Arguments must be different spinor fields\n");
    printf("Program aborted\n");
    exit(1);
  }

# if defined MPI
  xchange_spinorfield(Q);
# endif
   
  rho1.re=1.;
  rho1.im=g_mu;
  rho2.re=1.;
  rho2.im=-g_mu;

  /************************ loop over all lattice sites *************************/

  for (ix=0;ix<VOLUME;ix++){
    rr  = (spinor *) P +ix;
    s  = (spinor *) Q +ix;

    _complex_times_vector(r.s0,rho1,(*s).s0);
    _complex_times_vector(r.s1,rho1,(*s).s1);
    _complex_times_vector(r.s2,rho2,(*s).s2);
    _complex_times_vector(r.s3,rho2,(*s).s3);

    /******************************* direction +0 *********************************/

    iy=g_iup[ix][0];

    sp = (spinor *) Q +iy;
    up=&g_gauge_field[ix][0];
      
    _vector_add(psi,(*sp).s0,(*sp).s2);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_0, chi);

    _vector_add_assign(r.s0,psi);
    _vector_add_assign(r.s2,psi);

    _vector_add(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_0, chi);
    _vector_add_assign(r.s1,psi);
    _vector_add_assign(r.s3,psi);

    /******************************* direction -0 *********************************/

    iy=g_idn[ix][0];

    sm  = (spinor *) Q +iy;
    um=&g_gauge_field[iy][0];
      
    _vector_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_0, chi);
    _vector_add_assign(r.s0,psi);
    _vector_sub_assign(r.s2,psi);

    _vector_sub(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi,phase_0,chi);
    _vector_add_assign(r.s1,psi);
    _vector_sub_assign(r.s3,psi);

    /******************************* direction +1 *********************************/

    iy=g_iup[ix][1];

    sp = (spinor *) Q +iy;
    up=&g_gauge_field[ix][1];      
      
    _vector_i_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_1, chi);
    _vector_add_assign(r.s0,psi);
    _vector_i_sub_assign(r.s3,psi);

    _vector_i_add(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_1, chi);
    _vector_add_assign(r.s1,psi);
    _vector_i_sub_assign(r.s2,psi);

    /******************************* direction -1 *********************************/

    iy=g_idn[ix][1];

    sm = (spinor *) Q +iy;
    um=&g_gauge_field[iy][1];
      
    _vector_i_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_1, chi);
    _vector_add_assign(r.s0,psi);
    _vector_i_add_assign(r.s3,psi);

    _vector_i_sub(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_1, chi);
    _vector_add_assign(r.s1,psi);
    _vector_i_add_assign(r.s2,psi);

    /******************************* direction +2 *********************************/

    iy=g_iup[ix][2];

    sp = (spinor *) Q +iy;
    up=&g_gauge_field[ix][2];
      
    _vector_add(psi,(*sp).s0,(*sp).s3);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_2, chi);
    _vector_add_assign(r.s0,psi);
    _vector_add_assign(r.s3,psi);

    _vector_sub(psi,(*sp).s1,(*sp).s2);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_2, chi);
    _vector_add_assign(r.s1,psi);
    _vector_sub_assign(r.s2,psi);

    /******************************* direction -2 *********************************/

    iy=g_idn[ix][2];

    sm = (spinor *) Q +iy;
    um=&g_gauge_field[iy][2];
      
    _vector_sub(psi,(*sm).s0,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_2, chi);
    _vector_add_assign(r.s0,psi);
    _vector_sub_assign(r.s3,psi);

    _vector_add(psi,(*sm).s1,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_2, chi);
    _vector_add_assign(r.s1,psi);
    _vector_add_assign(r.s2,psi);

    /******************************* direction +3 *********************************/

    iy=g_iup[ix][3];

    sp = (spinor *) Q +iy;
    up=&g_gauge_field[ix][3];
      
    _vector_i_add(psi,(*sp).s0,(*sp).s2);
      
    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_3, chi);
    _vector_add_assign(r.s0,psi);
    _vector_i_sub_assign(r.s2,psi);

    _vector_i_sub(psi,(*sp).s1,(*sp).s3);

    _su3_multiply(chi,(*up),psi);

    _complex_times_vector(psi, phase_3, chi);
    _vector_add_assign(r.s1, psi);
    _vector_i_add_assign(r.s3, psi);

    /******************************* direction -3 *********************************/

    iy=g_idn[ix][3];

    sm = (spinor *) Q +iy;
    um=&g_gauge_field[iy][3];
      
    _vector_i_sub(psi,(*sm).s0,(*sm).s2);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_3, chi);
    _vector_add((*rr).s0, r.s0, psi);
    _vector_i_add((*rr).s2, r.s2, psi);

    _vector_i_add(psi,(*sm).s1,(*sm).s3);

    _su3_inverse_multiply(chi,(*um),psi);

    _complexcjg_times_vector(psi, phase_3, chi);
    _vector_add((*rr).s1, r.s1, psi);
    _vector_i_sub((*rr).s3, r.s3, psi);

    /******************************** end of loop *********************************/

  }
}

#endif

static char const rcsid[] = "$Id$";

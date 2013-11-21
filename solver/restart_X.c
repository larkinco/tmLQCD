/***********************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim based on the original code written
 * by Andreas Stathopoulos and Kostas Orginos and uses functions written in 
 * tmLQCD by Carsten Urbach 
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
 *21/11/2013:
 *In this version, we keep the BLAS implementation (it should be optimal) but
 *use multi-threading in case. In this case we use esize = (2*nev+1)*size of a long
 * spinor. This means that we have enough space to compute the result X*hVecs in one
 *go and store the result in rwork. This parameter is hard coded in incr_eigcg.c. 
 *We now add multi-threading by dividing the rows of X over the available threads.
 *
 ***********************************************************************/
/*******************************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 ******************************************************************************/
/**********************************************************
   <-------basisSize------>  <---restartSize-->
   |                      |  |                |
   |                      |  |                |
   |                      |  |                |
   |                      |  |   hVecs        |  
   |                      |  |                |
   |                      |  |                |
   |                      |  |                |
   |                      |  -----------------
   |       X              |
   |nLocal                |
   |                      |
   |ldx                   |
   <---------------------->
***************************************************/


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
#ifdef OMP
# include <omp.h>
#endif
#include "global.h"
#include "linalg_eo.h"
#include "su3.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"

#include "restart_X.h"
#define min(a, b) (a < b ? a : b)


void Zrestart_X(_Complex double  *X, int ldx, _Complex double  *hVecs, int nLocal, 
               int basisSize, int restartSize, _Complex double  *rwork, int rworkSize)
{
   /*old implementation
     ------------------
   char cN = 'N';
   int ONE = 1;
   int i, k;  // Loop variables 
   int AvailRows = min(rworkSize/restartSize, nLocal);
   _Complex double  tpone,tzero;
   tpone= +1.0e+00;  tzero=+0.0e+00;
   
   i = 0;

   while (i < nLocal) {
      // Block matrix multiply 
      _FT(zgemm)(&cN, &cN, &AvailRows, &restartSize, &basisSize, &tpone,
         &X[i], &ldx, hVecs, &basisSize, &tzero, rwork, &AvailRows ,1,1);

      // Copy the result in the desired location of X *
      for (k=0; k < restartSize; k++) {
         _FT(zcopy)(&AvailRows, &rwork[AvailRows*k],&ONE, &X[i+ldx*k],&ONE);
      }

      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }
   */
   /*
   char cN = 'N';
   int ONE = 1;
   int i, k;  // Loop variables 
   int AvailRows = min(rworkSize/restartSize, nLocal);
   _Complex double  tpone,tzero;
   tpone= +1.0e+00;  tzero=+0.0e+00;
   
   i = 0;

   while (i < nLocal) {
      // Block matrix multiply 
      _FT(zgemm)(&cN, &cN, &AvailRows, &restartSize, &basisSize, &tpone,
         &X[i], &ldx, hVecs, &basisSize, &tzero, rwork, &AvailRows ,1,1);

      // Copy the result in the desired location of X *
      for (k=0; k < restartSize; k++) {
         _FT(zcopy)(&AvailRows, &rwork[AvailRows*k],&ONE, &X[i+ldx*k],&ONE);
      }

      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }
   */

  
   
   char cN = 'N';
   int ONE = 1;
   int i, k;  // Loop variables 
   int AvailRows = min(rworkSize/restartSize, nLocal);
   _Complex double  tpone,tzero;
   tpone= +1.0e+00;  tzero=+0.0e+00;
 
   int nthreads=1;
   #ifdef OMP
      nthreads=omp_num_threads;
   #endif

   int nsize[nthreads]; //sizes of each partition

   int n1=nLocal/nthreads;

   int n2=nLocal-n1*nthreads;  //left over if nLocal is not integer multiple of nthreads;

   for(i=0; i<nthreads; i++)
      nsize[i]=n1;



   #ifdef OMP
     if(n2 >0)
     {
        //divide the left over points over the first n2 threads
        for(i=0; i<n2; i++)
           nsize[i] += 1;
     }    
   #endif

   int istart[nthreads]; //index of the first element for each partition
   istart[0]=0;
   for(i=1; i<nthreads; i++)
      istart[i]= istart[i-1]+nsize[i-1];



   #ifdef OMP
   #pragma omp for  
   for(i=0; i<nthreads; i++) {
      // Block matrix multiply 
      _FT(zgemm)(&cN, &cN, &nsize[i], &restartSize, &basisSize, &tpone,
         &X[istart[i]], &ldx, hVecs, &basisSize, &tzero, &rwork[istart[i]], &ldx,1,1);

      // Copy the result in the desired location of X *
      for (k=0; k < restartSize; k++) {
         _FT(zcopy)(&nsize[i], &rwork[istart[i]+ldx*k],&ONE, &X[istart[i]+ldx*k],&ONE);
      }
   }
   #endif

   

}




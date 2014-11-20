/*************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
 * Header file to interface with ARPACK package for solving large 
 * eigenvalue problems
 *
 * Author: Abdou M. Abdel-Rehim (amabdelrehim@gmail.com) 2014
 *************************************************************************/

#ifndef _ARPACK_H
#define _ARPACK_H

#include <complex.h>
#include "linalg/fortran.h"
#ifdef MPI
# include <mpi.h>
#endif




//ARPACK driver routines for computing eigenvectors (serial only) 
extern void _FT(znaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         _Complex double *resid, int *ncv, _Complex double *v, int *ldv, 
                         int *iparam, int *ipntr, _Complex double *workd, _Complex double *workl, 
                         int *lworkl, double *rwork, int *info );

extern void _FT(zneupd) (int *comp_evecs, char *cA, int *select, _Complex double *evals, 
                         _Complex double *v, int *ldv, _Complex double *sigma, _Complex double *workev, 
                         char *bmat, int *n, char *which, int *nev, double *tol, _Complex double *resid, 
                         int *ncv, _Complex double *v1, int *ldv1, int *iparam, int *ipntr, 
                         _Complex double *workd, _Complex double *workl, int *lworkl, double *rwork, int *info);

//PARPACK routines (parallel)
#ifdef MPI
extern void _FT(pznaupd) (MPI_Comm *comm, int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         _Complex double *resid, int *ncv, _Complex double *v, int *ldv, 
                         int *iparam, int *ipntr, _Complex double *workd, _Complex double *workl, 
                         int *lworkl, double *rwork, int *info );

extern void _FT(pzneupd) (MPI_Comm *comm, int *comp_evecs, char *cA, int *select, _Complex double *evals, 
                         _Complex double *v, int *ldv, _Complex double *sigma, _Complex double *workev, 
                         char *bmat, int *n, char *which, int *nev, double *tol, _Complex double *resid, 
                         int *ncv, _Complex double *v1, int *ldv1, int *iparam, int *ipntr, 
                         _Complex double *workd, _Complex double *workl, int *lworkl, double *rwork, int *info);
#endif

#endif

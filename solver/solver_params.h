/***************************************************************************
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
 ****************************************************************************/

/*****************************************************************************
 * Struct for passing the parameters for the solver. This should replace all 
 * solver related parameters and eliminate the need for global parameters as 
 * this struct will be a member of the operator struct.
 * 
 * A. M. Abdel-Rehim (a.abdel-rehim@cyi.ac.cy)
 * March, 17th, 2013
 ****************************************************************************/


#ifndef _SOLVER_PARAMS_H
#define _SOLVER_PARAMS_H

typedef struct {

  /********************************
   * Incremental EigCG parameters
   ********************************/

  int eigcg_nrhs;          /*total number of right-hand sides to be solved*/
  int eigcg_nrhs1;         /*The number of right-hand sides where we solve to tolerance tolsq1 
                             remaining systems will be solved to tolsq*/
  int eigcg_nev;           /*number of eigenvalues computed from a single right-hand side */
  int eigcg_vmax;          /*size of the search subspace for eigcg*/
  int eigcg_ldh;           /*total number of eigenvectors that will be computed and used in deflation */
  double eigcg_tolsq1;     /*squared tolerance for the first n1 systems */
  double eigcg_tolsq;      /*squared tolerance for the rest of the linear systems*/
  double eigcg_restolsq;   /*tolerance squared for restarting eigcg after eigenvectors has been computed
                             Typically this is the square root of the tolerance squared requested for the linear system.
                             Example, to solve the linear systems to squared residual 1e-16, one chooses eigcg_restolsq=1e-8 or smaller 
                             This will specify how many times deflated CG restaretd in the second phase (after eigenvectors has been computed)*/
  int eigcg_rand_guess_opt; /*set to 0 to use 0 initial guesses or non-zero values if you want to use random initial guess as a volume source */


  /**********************************
   * arpack_cg parameters
   **********************************/

   int    arpackcg_nrhs;          /*Number of right-hand sides to be solved*/ 
   int    arpackcg_nrhs1;         /*First number of right-hand sides to be solved using tolerance eps_sq1*/ 
   double arpackcg_eps_sq1;       /*Squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
   double arpackcg_eps_sq;        /*Squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
   double arpackcg_res_eps_sq;    /*Suqared tolerance for restarting cg */
   int    arpackcg_nev;           /*Number of eigenvectors to be computed by arpack*/
   int    arpackcg_ncv;           /*Size of the subspace used by arpack with the condition (nev+1) =< ncv =< 2*nev */
   double arpackcg_eig_tol;       /*Tolerance for computing eigenvalues with arpack */
   int    arpackcg_eig_maxiter;   /*Maximum number of iterations to be used by arpack*/

  /*************************************
   *chebychev polynomila preconditioner
   *for CG related paprameters
   *************************************/

   int cheb_k; //order of the polynomial used is k+1 and the lowest value is k=-1 which correspond to T_0
   double op_evmin;  //some estimate of the lowest eigenvalue
   double op_evmax;  //some estimate of the largest eigenvalue
   
   
} solver_params_t;



#endif


 

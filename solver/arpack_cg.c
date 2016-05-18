/*****************************************************************************
 * Copyright (C) 2014 Abdou M. Abdel-Rehim
 *
 * Deflating CG using eigenvectors computed using ARPACK
 * eigenvectors used correspond to those with smallest magnitude
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
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif

#include "global.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver_field.h"
#include "solver/arpack_cg.h"
#include "solver/eigen_topological.h"
# include "io/spinor.h"
# include "read_input.h"

int arpack_cg(

		//solver params
		const int N,             /*(IN) Number of lattice sites for this process*/
		const int nrhs,          /*(IN) Number of right-hand sides to be solved*/ 
		const int nrhs1,         /*(IN) First number of right-hand sides to be solved using tolerance eps_sq1*/ 
		spinor * const x,        /*(IN/OUT) initial guess on input, solution on output for this RHS*/
		spinor * const b,        /*(IN) right-hand side*/
		matrix_mult f,           /*(IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply*/
		const double eps_sq1,    /*(IN) squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
		const double eps_sq,     /*(IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
		const double res_eps_sq, /*(IN) suqared tolerance for restarting cg */
		const int rel_prec,      /*(IN) 0 for using absoute error for convergence
					   1 for using relative error for convergence*/
		const int maxit,         /*(IN) Maximum allowed number of iterations to solution for the linear system*/

		//parameters for arpack
		const int nev,                 /*(IN) number of eigenvectors to be computed by arpack*/
		const int ncv,                 /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
		double arpack_eig_tol,         /*(IN) tolerance for computing eigenvalues with arpack */
		int arpack_eig_maxiter,        /*(IN) maximum number of iterations to be used by arpack*/
		int kind,                      /*(IN) 0 for eigenvalues with smallest real part "SR"
						 1 for eigenvalues with largest real part "LR"
						 2 for eigenvalues with smallest absolute value "SM"
						 3 for eigenvalues with largest absolute value "LM"
						 4 for eigenvalues with smallest imaginary part "SI"
						 5 for eigenvalues with largest imaginary part  "LI"*/
		int comp_evecs,                /*(IN) 0 don't compute the eiegnvalues and their residuals of the original system 
						 1 compute the eigenvalues and the residuals for the original system (the orthonormal baiss
						 still be used in deflation and they are not overwritten).*/
		int acc,                       /*(IN) 0 no polynomial acceleration
						 1 use polynomial acceleration*/
		int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
		double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
		double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
		int read_basis,               /*(IN) 0 compute deflation basis using arpack, 1 read deflation basis from disk */
		int store_basis,              /*(IN) option to store basis vectors to disk such that they can be read later
						0 don't store basis vectors
						1 store basis vectors */
		char *basis_fname,            /*(IN)file name used to read/store the basis vectors
						file names will be of the format
						basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
		int basis_prec,               /*(IN)precision used to write the basis vectors
						0 single precision
						1 double precision*/
		char *arpack_logfile,           /*(IN) file name to be used for printing out debugging information from arpack*/		
		int top_suscept_exact_flag				
)
{ 

	//Static variables and arrays.
	static int ncurRHS=0;                  /* current number of the system being solved */                   
	static void *_ax,*_r,*_tmps1,*_tmps2,*_vec_even,*_vec_odd;                  
	static spinor *ax,*r,*tmps1,*tmps2,*vec_even,*vec_odd;                  
	static _Complex double *evecs,*evals,*H,*HU,*Hinv,*initwork,*tmpv1;
	static _Complex double *zheev_work;
	static double *hevals,*zheev_rwork;
	_Complex double* ef_evs=NULL;
	static int *IPIV; 
	static int info_arpack=0;
	static int nconv=0; //number of converged eigenvectors as returned by arpack
	int i,j,tmpsize;
	char cV='V',cN='N', cU='U';   
	int ONE=1;
	int zheev_lwork,zheev_info;
	_Complex double c1,c2,tpone=1.0,tzero=0.0;
	double d1,d2,d3;
	double et1,et2;  //timing variables

	int prec,status,rstat;

	int parallel;        /* for parallel processing of the scalar products */
#ifdef MPI
	parallel=1;
#else
	parallel=0;
#endif

	/* leading dimension for spinor vectors */
	int LDN;
	if(N==VOLUME)
		LDN = VOLUMEPLUSRAND;
	else
		LDN = VOLUMEPLUSRAND/2; 

	/*-------------------------------------------------------------
	//if this is the first right hand side, allocate memory, 
	//call arpack, and compute resiudals of eigenvectors if needed
	//-------------------------------------------------------------*/ 
	if(ncurRHS==0){ 
#if (defined SSE || defined SSE2 || defined SSE3)
		_ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_ax==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _ax inside arpack_cg.\n");
			exit(1);
		}
		else
		{ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

		_r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_r==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _r inside arpack_cg.\n");
			exit(1);
		}
		else
		{r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

		_tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_tmps1==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _tmps1 inside arpack_cg.\n");
			exit(1);
		}
		else
		{tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

		_tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_tmps2==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _tmps2 inside arpack_cg.\n");
			exit(1);
		}
		else
		{tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


		_vec_even = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_vec_even==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _vec_even inside arpack_cg.\n");
			exit(1);
		}
		else
		{vec_even  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


		_vec_odd = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
		if(_vec_odd==NULL)
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for _vec_odd inside arpack_cg.\n");
			exit(1);
		}
		else
		{vec_odd  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}

#else
		ax = (spinor *) malloc(LDN*sizeof(spinor));
		r  = (spinor *) malloc(LDN*sizeof(spinor));
		tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
		tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
		vec_even = (spinor *) malloc(LDN*sizeof(spinor));  
		vec_odd  = (spinor *) malloc(LDN*sizeof(spinor));

		if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) || (vec_even==NULL) || (vec_odd==NULL) )
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for ax,r,tmps1,tmps2,vec_even,vec_odd inside arpack_cg.\n");
			exit(1);
		}
#endif


		evecs = (_Complex double *) malloc(ncv*12*N*sizeof(_Complex double)); //note: no extra buffer 
		evals = (_Complex double *) malloc(ncv*sizeof(_Complex double)); 
		tmpv1 = (_Complex double *) malloc(12*N*sizeof(_Complex double));

		if((evecs == NULL)  || (evals==NULL) || (tmpv1==NULL))
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for evecs and evals inside arpack_cg.\n");
			exit(1);
		}

		if(read_basis)
		{
			et1=gettime();
			nconv=nev;
			for(j=0; j<nev; j++)
			{
				char filename[500],*header_type=NULL; 
				READER *reader=NULL;
				uint64_t bytes;
				//sprintf(filename, "ev.%04d.%05d", nstore, i);
				sprintf(filename, "%s.%05d", basis_fname, j);
				construct_reader(&reader,filename); 
				DML_Checksum checksum;

				/* Find the desired binary data*/
				while ((status = ReaderNextRecord(reader)) != LIME_EOF) {
					if (status != LIME_SUCCESS){
						fprintf(stderr, "ReaderNextRecord returned status %d.\n", status);
						break;
					}
					header_type = ReaderType(reader);
					if (strcmp("scidac-binary-data", header_type) == 0) {
						break;
					}
				}

				if (status == LIME_EOF) {
					fprintf(stderr, "Unable to find requested LIME record scidac-binary-data in file %s.\nEnd of file reached before record was found.\n", filename);
					return(-5);
				}

				bytes = ReaderBytes(reader);

				if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor)) {
					prec = 64;
				}
				else {
					if ((int)bytes == LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * sizeof(spinor) / 2) {
						prec = 32;
					}
					else {
						fprintf(stderr, "Length of scidac-binary-data record in %s does not match input parameters.\n", filename);
						fprintf(stderr, "Found %d bytes.\n", bytes);
						return(-6);
					}
				}

				if (g_cart_id == 0 && g_debug_level >= 0) {
					printf("# %s precision read (%d bits).\n", (prec == 64 ? "Double" : "Single") ,prec);
				}

				if( (rstat = read_binary_spinor_data(vec_even,vec_odd, reader, &checksum)) != 0) {
					fprintf(stderr, "read_binary_spinor_data failed with return value %d", rstat);
					return(-7);
				}

				if(N==VOLUME) //solving the full system
				{
					convert_eo_to_lexic(r,vec_even,vec_odd);
					assign_spinor_to_complex(&evecs[j*12*N],r,N);
				}else{  //solving the eo preconditioned systems (only odd part of the eigenvector is computed)
					assign_spinor_to_complex(&evecs[j*12*N],vec_odd,N);
				}


				if (g_cart_id == 0 && g_debug_level >= 0) {
					printf("# Scidac checksums for DiracFermion field %s:\n", filename);
					printf("#   Calculated            : A = %#x B = %#x.\n", checksum.suma, checksum.sumb);
					printf("# No Scidac checksum was read from headers, unable to check integrity of file.\n");
				}

				destruct_reader(reader);
			} //for(j=0;...)
			et2=gettime();
			if(g_proc_id == g_stdio_proc)
				fprintf(stdout,"Finished reading deflation basis in %e seconds\n",et2-et1); 
		}
		else
		{
			et1=gettime();
			evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv,arpack_logfile);
			et2=gettime();


			if(info_arpack != 0){ //arpack didn't converge
				if(g_proc_id == g_stdio_proc)
					fprintf(stderr,"WARNING: ARPACK didn't converge. exiting..\n");
				return -1;
			}

			if(g_proc_id == g_stdio_proc)
			{
				fprintf(stdout,"ARPACK has computed %d eigenvectors\n",nconv);
				fprintf(stdout,"ARPACK time: %+e\n",et2-et1);
			}
		}

		H        = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
		HU       = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
		Hinv     = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
		initwork = (_Complex double *) malloc(nconv*sizeof(_Complex double)); 
		IPIV     = (int *) malloc(nconv*sizeof(int));
		zheev_lwork = 3*nconv;
		zheev_work  = (_Complex double *) malloc(zheev_lwork*sizeof(_Complex double));
		zheev_rwork = (double *) malloc(3*nconv*sizeof(double));
		hevals      = (double *) malloc(nconv*sizeof(double));

		if((H==NULL) || (HU==NULL) || (Hinv==NULL) || (initwork==NULL) || (IPIV==NULL) || (zheev_lwork==NULL) || (zheev_rwork==NULL) || (hevals==NULL))
		{
			if(g_proc_id == g_stdio_proc)
				fprintf(stderr,"insufficient memory for H, HU, Hinv, initwork, IPIV, zheev_lwork, zheev_rwork, hevals inside arpack_cg.\n");
			exit(1);
		}

		//compute the elements of the hermitian matrix H 
		//leading dimension is nconv and active dimension is nconv
		for(i=0; i<nconv; i++)
		{
			assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
			f(ax,r);
			c1 = scalar_prod(r,ax,N,parallel);
			H[i+nconv*i] = creal(c1);  //diagonal should be real
			for(j=i+1; j<nconv; j++)
			{
				assign_complex_to_spinor(r,&evecs[j*12*N],12*N);
				c1 = scalar_prod(r,ax,N,parallel);
				H[j+nconv*i] = c1;
				H[i+nconv*j] = conj(c1); //enforce hermiticity
			}
		}

		//compute Ritz values and Ritz vectors if needed
		if( (nconv>0) && (comp_evecs !=0))
		{
			/* copy H into HU */
			tmpsize=nconv*nconv;
			_FT(zcopy)(&tmpsize,H,&ONE,HU,&ONE);

			/* compute eigenvalues and eigenvectors of HU*/
			//SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,INFO )
			_FT(zheev)(&cV,&cU,&nconv,HU,&nconv,hevals,zheev_work,&zheev_lwork,zheev_rwork,&zheev_info,1,1); //hevals found at this point

			if(zheev_info != 0)
			{
				if(g_proc_id == g_stdio_proc) 
				{
					fprintf(stderr,"Error in ZHEEV:, info =  %d\n",zheev_info); 
					fflush(stderr);
				}
				exit(1);
			}

			//If you want to replace the schur (orthonormal) basis by eigen basis
			//use something like this. It is better to use the schur basis because
			//they are better conditioned. Use this part only to get the eigenvalues
			//and their resduals for the operator (D^\daggerD)
			//esize=(ncv-nconv)*12*N;
			//Zrestart_X(evecs,12*N,HU,12*N,nconv,nconv,&evecs[nconv*N],esize);

			/* compute residuals and print out results */

			

//			printf("The value of top flag is %d",top_suscept_exact_flag);
			//printf("WE ARE HERE");
			if(N==VOLUME&&top_suscept_exact_flag==1)
			{
				ef_evs =  malloc(nconv*12*N*sizeof(_Complex double)); //note: no extra buffer 
				//CHECK NOT NULL

				if(ef_evs==NULL)
				{
					if(g_proc_id == g_stdio_proc)
						fprintf(stderr,"insufficient memory for ef_evs inside invert.\n");
					exit(1);
				}
			}


			if(g_proc_id == g_stdio_proc)
			{fprintf(stdout,"Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
				fprintf(stdout,"=============================================================\n");
				fflush(stdout);}

			for(i=0; i<nconv; i++)
			{
				tmpsize=12*N;
				_FT(zgemv)(&cN,&tmpsize,&nconv,&tpone,evecs,&tmpsize,
						&HU[i*nconv],&ONE,&tzero,tmpv1,&ONE,1);


				assign_complex_to_spinor(r,tmpv1,12*N);

				
				if(N==VOLUME&&top_suscept_exact_flag==1)
				{
					memcpy(ef_evs+tmpsize*i,tmpv1,tmpsize*sizeof(_Complex double)); //memcopying all individual eigenvectors into our full converged eigenvector array used in topological calculation. This is done to require minimum changing of the working code to minimise testing time. 
				}


				if(store_basis){

					WRITER *writer = NULL;
					int append = 0;
					char fname[256];	
					paramsPropagatorFormat *format = NULL;
					int precision;
					int numb_flavs = 1;


					if(basis_prec==0)
						precision = 32;
					else
						precision = 64;

					sprintf(fname, "%s.%05d", basis_fname, i);
					construct_writer(&writer, fname, append);

					format = construct_paramsPropagatorFormat(precision, numb_flavs);
					write_propagator_format(writer, format);
					free(format);	    

					//memset(tmpv2, '\0', size*sizeof(_Complex double));
					//assign_complex_to_spinor(tmps2,tmpv2,size);
					int status;
					if(N==VOLUME){ //solving the full system
						convert_lexic_to_eo(vec_even,vec_odd,r);
						status = write_spinor(writer, &vec_even, &vec_odd, numb_flavs, precision);
					}else{
						zero_spinor_field(vec_even,N);
						status = write_spinor(writer, &vec_even, &r, numb_flavs, precision);
					}
					destruct_writer(writer);

				} //if(store_basis)...

				d1=square_norm(r,N,parallel);

				f(ax,r);

				mul_r(tmps1,hevals[i],r,N);

				diff(tmps2,ax,tmps1,N);

				d2= square_norm(tmps2,N,parallel);

				d3= sqrt(d2/d1);

				if(g_proc_id == g_stdio_proc)
				{fprintf(stdout,"Eval[%06d]: %22.15E rnorm: %22.15E\n", i, hevals[i], d3); fflush(stdout);}
			}

			/////////////////////////////////////////////// 
			/////////////////////////////////////////////// 
			//CHECK EVEN ODD AND TOPOLOGICA SUSCEPT FLAG AND IF BOTH TRUE RUN
			/*
			//char top_filename[] = "/homec/ecy00/ecy002/Conor/Runs/top_output.txt";
			int top_suscept_exact_flag=1;//temp for testing
			//printf("WE ARE HERE");
			if(N==VOLUME&&top_suscept_exact_flag==1)
			{

			eigen_field eigf;//CHECK EIGN NUM
			eigen_field_init(evecs,hevals,nconv,N,&eigf);

			//	_Complex double final_answer[400];
			int n_eig=400;
			_Complex double* top_suscept= malloc(sizeof(_Complex double)*n_eig);
			if(top_suscept==NULL)
			{
			if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _ax inside arpack.\n");
			exit(1);
			}

			for(int i =1;i<n_eig-20;i++)
			{
			top_suscept[i-1] = top_suscept_subset(&eigf,i);

			double x=0.005501;
			top_suscept[i-1] = top_suscept[i-1]*pow(x,6);
			if(g_proc_id == g_stdio_proc)
			{ 
			printf("WE ARE INSIDE");	
			printf("%d %.10e %.10e \n",i,creal(top_suscept[i-1]),cimag(top_suscept[i-1]));

			}		
			}
			}
			*/


		//	printf("THE value of top flag is %d.",top_suscept_exact_flag);
			if(N==VOLUME&&top_suscept_exact_flag==1)
			{

				eigen_field eigf;//CHECK EIGN NUM
				eigen_field_init(ef_evs,hevals,nconv,N,&eigf);//Change evecs to new eigenvector array

				//	_Complex double final_answer[400];
				//	int n_eig=nconv;
				if(g_proc_id == g_stdio_proc)
				{ 
					printf("WE ARE INSide top if statement, before allocation and calculation");	
				}
				_Complex double* top_suscept= malloc(sizeof(_Complex double)*nconv);
				if(top_suscept==NULL)
				{
					if(g_proc_id == g_stdio_proc)
						fprintf(stderr,"insufficient memory for top_suscept array.\n");
					exit(1);
				}

				for(int i =1;i<=nconv;i++)
				{
					top_suscept[i-1] = top_suscept_subset(&eigf,i);

					double x=0.005501;
					top_suscept[i-1] = top_suscept[i-1]*pow(x,6);
					if(g_proc_id == g_stdio_proc)
					{ 
						//										printf("WE ARE INSIDE");	
						printf("%d %.10e %.10e \n",i,creal(top_suscept[i-1]),cimag(top_suscept[i-1]));

					}		
				_Complex double temp= top_suscept_product_print(&eigf,nconv);	
				}
				free(top_suscept);
			}


			/////////////////////////////////////////////// 
			/////////////////////////////////////////////// 
		}//if( (nconv_arpack>0) && (comp_evecs !=0))
	} //if(ncurRHS==0)

	double eps_sq_used,restart_eps_sq_used;  //tolerance squared for the linear system

	double cur_res; //current residual squared

	/*increment the RHS counter*/
	ncurRHS = ncurRHS +1; 

	//set the tolerance to be used for this right-hand side 
	if(ncurRHS > nrhs1){
		eps_sq_used = eps_sq;
	}
	else{
		eps_sq_used = eps_sq1;
	}

	if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
		fprintf(stdout, "System %d, eps_sq %e\n",ncurRHS,eps_sq_used); 
		fflush(stdout);
	} 

	/*---------------------------------------------------------------*/
	/* Call init-CG until this right-hand side converges             */
	/*---------------------------------------------------------------*/
	double wt1,wt2,wE,wI;
	double normsq,tol_sq;
	int flag,maxit_remain,numIts,its;
	int info_lapack;

	wE = 0.0; wI = 0.0;     /* Start accumulator timers */
	flag = -1;    	  /* System has not converged yet */
	maxit_remain = maxit;   /* Initialize Max and current # of iters   */
	numIts = 0;  
	restart_eps_sq_used=res_eps_sq;

	while( flag == -1 )
	{

		if(nconv > 0)
		{
			/* --------------------------------------------------------- */
			/* Perform init-CG with evecs vectors                        */
			/* xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		   */
			/* --------------------------------------------------------- */
			wt1 = gettime();

			/*r0=b-Ax0*/
			f(ax,x); /*ax = A*x */
			diff(r,b,ax,N);  /* r=b-A*x */

			/* x = x + evecs*inv(H)*evecs'*r */
			for(int i=0; i < nconv; i++)
			{
				assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
				initwork[i]= scalar_prod(tmps1,r,N,parallel);
			}

			/* solve the linear system H y = c */
			tmpsize=nconv*nconv;
			_FT(zcopy) (&tmpsize,H,&ONE,Hinv,&ONE); /* copy H into Hinv */
			//SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
			_FT(zgesv) (&nconv,&ONE,Hinv,&nconv,IPIV,initwork,&nconv,&info_lapack);

			if(info_lapack != 0)
			{
				if(g_proc_id == g_stdio_proc) {
					fprintf(stderr, "Error in ZGESV:, info =  %d\n",info_lapack); 
					fflush(stderr);
				}
				exit(1);
			}

			/* x = x + evecs*inv(H)*evecs'*r */
			for(i=0; i<nconv; i++)
			{
				assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
				assign_add_mul(x,tmps1,initwork[i],N);
			}

			/* compute elapsed time and add to accumulator */

			wt2 = gettime();
			wI = wI + wt2-wt1;

		}/* if(nconv > 0) */


		//which tolerance to use
		if(eps_sq_used > restart_eps_sq_used)
		{
			tol_sq = eps_sq_used;
			flag   = 1; //shouldn't restart again
		}
		else
		{
			tol_sq = restart_eps_sq_used;
		}

		wt1 = gettime();
		its = cg_her(x,b,maxit_remain,tol_sq,rel_prec,N,f); 

		wt2 = gettime();

		wE = wE + wt2-wt1;

		//check convergence
		if(its == -1)
		{
			//cg didn't converge
			if(g_proc_id == g_stdio_proc) {
				fprintf(stderr, "CG didn't converge within the maximum number of iterations in arpack_cg. Exiting...\n"); 
				fflush(stderr);
				exit(1);

			}
		} 
		else
		{
			numIts += its;   
			maxit_remain = maxit - numIts; //remaining number of iterations
			restart_eps_sq_used = restart_eps_sq_used*res_eps_sq; //prepare for the next restart
		}

	}
	/* end while (flag ==-1)               */

	/* ---------- */
	/* Reporting  */
	/* ---------- */
	/* compute the exact residual */
	f(ax,x); /* ax= A*x */
	diff(r,b,ax,N);  /* r=b-A*x */	
	normsq=square_norm(r,N,parallel);
	if(g_debug_level > 0 && g_proc_id == g_stdio_proc)
	{
		fprintf(stdout, "For this rhs:\n");
		fprintf(stdout, "Total initCG Wallclock : %+e\n", wI);
		fprintf(stdout, "Total cg Wallclock : %+e\n", wE);
		fprintf(stdout, "Iterations: %-d\n", numIts); 
		fprintf(stdout, "Actual Resid of LinSys  : %+e\n",normsq);
	}


	//free memory if this was your last system to solve
	if(ncurRHS == nrhs){
#if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
		free(_ax);  free(_r);  free(_tmps1); free(_tmps2); free(_vec_even); free(_vec_odd);
#else
		free(ax); free(r); free(tmps1); free(tmps2); free(vec_even); free(vec_odd);
#endif
		free(evecs); free(evals); free(H); free(HU); free(Hinv);
		free(initwork); free(tmpv1); free(zheev_work);
		free(hevals); free(zheev_rwork); free(IPIV);free(ef_evs);
	}


	return numIts;
}


int arpack(

		//solver params
		const int N,             /*(IN) Number of lattice sites for this process*/
		matrix_mult f,           /*(IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply*/
		//parameters for arpack
		const int nev,                 /*(IN) number of eigenvectors to be computed by arpack*/
		const int ncv,                 /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
		double arpack_eig_tol,         /*(IN) tolerance for computing eigenvalues with arpack */
		int arpack_eig_maxiter,        /*(IN) maximum number of iterations to be used by arpack*/
		int kind,                      /*(IN) 0 for eigenvalues with smallest real part "SR"
						 1 for eigenvalues with largest real part "LR"
						 2 for eigenvalues with smallest absolute value "SM"
						 3 for eigenvalues with largest absolute value "LM"
						 4 for eigenvalues with smallest imaginary part "SI"
						 5 for eigenvalues with largest imaginary part  "LI"*/
		int acc,                       /*(IN) 0 no polynomial acceleration
						 1 use polynomial acceleration*/
		int cheb_k,                    /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
		double emin,                      /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
		double emax,                      /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
		int store_basis,              /*(IN) option to store basis vectors to disk such that they can be read later
						0 don't store basis vectors
						1 store basis vectors */
		char *basis_fname,            /*(IN)file name used to read/store the basis vectors
						file names will be of the format
						basis_fname.xxxxx where xxxxx will be the basis vector number with leading zeros */
		int basis_prec,               /*(IN)precision used to write the basis vectors
						0 single precision
						1 double precision*/
		char *arpack_logfile,           /*(IN) file name to be used for printing out debugging information from arpack*/
		int top_suscept_exact_flag				
		)
{ 

	//Static variables and arrays.
	static void *_ax,*_r,*_tmps1,*_tmps2,*_vec_even,*_vec_odd;                  
	static spinor *ax,*r,*tmps1,*tmps2,*vec_even,*vec_odd;                  
	static _Complex double *evecs,*evals,*H,*HU,*Hinv,*initwork,*tmpv1;
	static _Complex double *zheev_work;
	static double *hevals,*zheev_rwork;
	static int *IPIV; 
	static int info_arpack=0;
	static int nconv=0; //number of converged eigenvectors as returned by arpack
	int i,j,tmpsize;
	char cV='V',cN='N', cU='U';   
	int ONE=1;
	int zheev_lwork,zheev_info;
	_Complex double c1,c2,tpone=1.0,tzero=0.0;
	double d1,d2,d3;
	double et1,et2;  //timing variables
	_Complex double* ef_evs=NULL;
	int prec,status,rstat;

	int parallel;        /* for parallel processing of the scalar products */
#ifdef MPI
	parallel=1;
#else
	parallel=0;
#endif

	/* leading dimension for spinor vectors */
	int LDN;
	if(N==VOLUME)
		LDN = VOLUMEPLUSRAND;
	else
		LDN = VOLUMEPLUSRAND/2; 

	/*-------------------------------------------------------------
	//allocate memory 
	//-------------------------------------------------------------*/  
#if (defined SSE || defined SSE2 || defined SSE3)
	_ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_ax==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _ax inside arpack.\n");
		exit(1);
	}
	else
	{ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

	_r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_r==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _r inside arpack.\n");
		exit(1);
	}
	else
	{r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

	_tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_tmps1==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _tmps1 inside arpack.\n");
		exit(1);
	}
	else
	{tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

	_tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_tmps2==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _tmps2 inside arpack.\n");
		exit(1);
	}
	else
	{tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


	_vec_even = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_vec_even==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _zero_spinor inside arpack.\n");
		exit(1);
	}
	else
	{vec_even  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}

	_vec_odd = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
	if(_vec_odd==NULL)
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for _zero_spinor inside arpack.\n");
		exit(1);
	}
	else
	{vec_odd  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}


#else
	ax = (spinor *) malloc(LDN*sizeof(spinor));
	r  = (spinor *) malloc(LDN*sizeof(spinor));
	tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
	tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
	vec_even = (spinor *) malloc(LDN*sizeof(spinor));
	vec_odd  = (spinor *) malloc(LDN*sizeof(spinor));

	if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) || (vec_even==NULL) || (vec_odd==NULL))
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for ax,r,tmps1,tmps2,zero_spinor inside arpack.\n");
		exit(1);
	}
#endif

	evecs = (_Complex double *) malloc(ncv*12*N*sizeof(_Complex double)); //note: no extra buffer 
	evals = (_Complex double *) malloc(ncv*sizeof(_Complex double)); 
	tmpv1 = (_Complex double *) malloc(12*N*sizeof(_Complex double));

	if((evecs == NULL)  || (evals==NULL) || (tmpv1==NULL))
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for evecs and evals inside arpack.\n");
		exit(1);
	}

	et1=gettime();
	evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv,arpack_logfile);
	et2=gettime();


	if(info_arpack != 0){ //arpack didn't converge
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"WARNING: ARPACK didn't converge. exiting..\n");
		exit(1);
	}

	if(g_proc_id == g_stdio_proc)
	{
		fprintf(stdout,"ARPACK has computed %d eigenvectors\n",nconv);
		fprintf(stdout,"ARPACK time: %+e\n",et2-et1);
	}

	H        = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
	HU       = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
	Hinv     = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
	initwork = (_Complex double *) malloc(nconv*sizeof(_Complex double)); 
	IPIV     = (int *) malloc(nconv*sizeof(int));
	zheev_lwork = 3*nconv;
	zheev_work  = (_Complex double *) malloc(zheev_lwork*sizeof(_Complex double));
	zheev_rwork = (double *) malloc(3*nconv*sizeof(double));
	hevals      = (double *) malloc(nconv*sizeof(double));

	if((H==NULL) || (HU==NULL) || (Hinv==NULL) || (initwork==NULL) || (IPIV==NULL) || (zheev_lwork==NULL) || (zheev_rwork==NULL) || (hevals==NULL))
	{
		if(g_proc_id == g_stdio_proc)
			fprintf(stderr,"insufficient memory for H, HU, Hinv, initwork, IPIV, zheev_lwork, zheev_rwork, hevals inside arpack.\n");
		exit(1);
	}

	//compute the elements of the hermitian matrix H 
	//leading dimension is nconv and active dimension is nconv
	for(i=0; i<nconv; i++)
	{
		assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
		f(ax,r);
		c1 = scalar_prod(r,ax,N,parallel);
		H[i+nconv*i] = creal(c1);  //diagonal should be real
		for(j=i+1; j<nconv; j++)
		{
			assign_complex_to_spinor(r,&evecs[j*12*N],12*N);
			c1 = scalar_prod(r,ax,N,parallel);
			H[j+nconv*i] = c1;
			H[i+nconv*j] = conj(c1); //enforce hermiticity
		}
	}

	//compute Ritz values and Ritz vectors if needed
	if( (nconv>0))
	{
		/* copy H into HU */
		tmpsize=nconv*nconv;
		_FT(zcopy)(&tmpsize,H,&ONE,HU,&ONE);

		/* compute eigenvalues and eigenvectors of HU*/
		//SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,INFO )
		_FT(zheev)(&cV,&cU,&nconv,HU,&nconv,hevals,zheev_work,&zheev_lwork,zheev_rwork,&zheev_info,1,1);

		if(zheev_info != 0)
		{
			if(g_proc_id == g_stdio_proc) 
			{
				fprintf(stderr,"Error in ZHEEV:, info =  %d\n",zheev_info); 
				fflush(stderr);
			}
			exit(1);
		}

		//If you want to replace the schur (orthonormal) basis by eigen basis
		//use something like this. It is better to use the schur basis because
		//they are better conditioned. Use this part only to get the eigenvalues
		//and their resduals for the operator (D^\daggerD)
		//esize=(ncv-nconv)*12*N;
		//Zrestart_X(evecs,12*N,HU,12*N,nconv,nconv,&evecs[nconv*N],esize);

		/* compute residuals and print out results */

		if(g_proc_id == g_stdio_proc)
		{fprintf(stdout,"Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
			fprintf(stdout,"=============================================================\n");
			fflush(stdout);}

		//If top_calc
		//	create a new eigenvector array of size nconv*12*N
			int top_suscept_exact_flag=1;//temp for testing
//		printf("The value of top flag is %d",top_suscept_exact_flag);
		//printf("WE ARE HERE");
		if(N==VOLUME&&top_suscept_exact_flag==1)
		{
			ef_evs =  malloc(nconv*12*N*sizeof(_Complex double)); //note: no extra buffer 
			//CHECK NOT NULL

			if(ef_evs==NULL)
			{
				if(g_proc_id == g_stdio_proc)
					fprintf(stderr,"insufficient memory for ef_evs inside invert.\n");
				exit(1);
			}
		}


		for(i=0; i<nconv; i++)
		{
			tmpsize=12*N;//size of an eigenvector
			_FT(zgemv)(&cN,&tmpsize,&nconv,&tpone,evecs,&tmpsize,
					&HU[i*nconv],&ONE,&tzero,tmpv1,&ONE,1); //tmpv1 IS the output, rest are input, multiplies, tempv1 = evecs*HU[i*nconv] + 0

			if(N==VOLUME&&top_suscept_exact_flag==1)
			{
				memcpy(ef_evs+tmpsize*i,tmpv1,tmpsize*sizeof(_Complex double)); //memcopying all individual eigenvectors into our full converged eigenvector array used in topological calculation. This is done to require minimum changing of the working code to minimise testing time. 
			}

			assign_complex_to_spinor(r,tmpv1,12*N);

			if(store_basis){
				WRITER *writer = NULL;
				int append = 0;
				char fname[256];	
				paramsPropagatorFormat *format = NULL;
				int precision;
				int numb_flavs = 1;


				if(basis_prec==0)
					precision = 32;
				else
					precision = 64;

				sprintf(fname, "%s.%05d", basis_fname, i);
				construct_writer(&writer, fname, append);

				format = construct_paramsPropagatorFormat(precision, numb_flavs);
				write_propagator_format(writer, format);
				free(format);	    

				//memset(tmpv2, '\0', size*sizeof(_Complex double));
				//assign_complex_to_spinor(tmps2,tmpv2,size);
				int status;
				if(N==VOLUME){  //full system
					convert_lexic_to_eo(vec_even,vec_odd,r);
					status = write_spinor(writer, &vec_even, &vec_odd, numb_flavs, precision);
				}
				else{ //EO case
					zero_spinor_field(vec_even,N);
					status = write_spinor(writer, &vec_even, &r, numb_flavs, precision);
				}
				destruct_writer(writer);

			} //if(store_basis)...

			d1=square_norm(r,N,parallel);

			f(ax,r);

			mul_r(tmps1,hevals[i],r,N);

			diff(tmps2,ax,tmps1,N);

			d2= square_norm(tmps2,N,parallel);

			d3= sqrt(d2/d1);

			if(g_proc_id == g_stdio_proc)
			{fprintf(stdout,"Eval[%06d]: %22.15E rnorm: %22.15E\n", i, hevals[i], d3); fflush(stdout);}
		}

		/////////////////////////////////////////////// 
		/////////////////////////////////////////////// 
		//CHECK EVEN ODD AND TOPOLOGICA SUSCEPT FLAG AND IF BOTH TRUE RUN
		/*char top_filename[] = "/homec/ecy00/ecy002/Conor/Runs/top_output.txt";
		  int top_suscept_exact_flag=1;//temp for testing
		//	printf("WE ARE HERE");
		if(N==VOLUME&&top_suscept_exact_flag==1)
		{
		//	printf("WE ARE INSIDE");

		if(g_proc_id == g_stdio_proc)
		{ 	
		//eigen_helloworld();
		printf("WE ARE INSIDE %d", N);
		FILE * outputfile;
		outputfile = fopen (top_filename, "r+");
		fprintf(outputfile,"helloworld");
		fclose(outputfile);
		}		
		}
		*/

		//	int top_suscept_exact_flag=1;//temp for testing
		//printf("WE ARE HERE");
		//printf("THE value of top flag is %d.",top_suscept_exact_flag);
		if(N==VOLUME&&top_suscept_exact_flag==1)
		{

			eigen_field eigf;//CHECK EIGN NUM
			eigen_field_init(ef_evs,hevals,nconv,N,&eigf);//Change evecs to new eigenvector array

			//	_Complex double final_answer[400];
			//	int n_eig=nconv;
			if(g_proc_id == g_stdio_proc)
			{ 
				printf("WE ARE INSide top if statement, before allocation and calculation");	
			}
			_Complex double* top_suscept= malloc(sizeof(_Complex double)*nconv);
			if(top_suscept==NULL)
			{
				if(g_proc_id == g_stdio_proc)
					fprintf(stderr,"insufficient memory for top_suscept array.\n");
				exit(1);
			}

			for(int i =1;i<=nconv;i++)
			{
				top_suscept[i-1] = top_suscept_subset(&eigf,i);

				double x=0.005501;
				top_suscept[i-1] = top_suscept[i-1]*pow(x,6);
				
			//	_Complex double Q_sqred2=Q_squared_2(&eigf, i-1);
			//	_Complex double Q_sqred1=Q_squared_1(&eigf, i-1);
				//_Complex double renrm_top_sus=renorm_top_suscept_alt(&eigf,i);
				if(g_proc_id == g_stdio_proc)
				{ 
					//									printf("WE ARE INSIDE");	
					printf("%d %.10e %.10e \n",i,creal(top_suscept[i-1]),cimag(top_suscept[i-1]));
			//		printf("n Q_sqred1 %d %.10e %.10e \n",i,creal(Q_sqred1),cimag(Q_sqred1));
			//		printf("n Q_sqred2 %d %.10e %.10e \n",i,creal(Q_sqred2),cimag(Q_sqred2));
				//	printf("n renrm_top_sus %d %.10e %.10e \n",i,creal(renrm_top_sus),cimag(renrm_top_sus));
				}		
			}
			_Complex double temp= top_suscept_product_print(&eigf,nconv);
			free(top_suscept);
		}
		/////////////////////////////////////////////// 
		/////////////////////////////////////////////// 


	}//if( (nconv_arpack>0) && (comp_evecs !=0))


	//free memory 
#if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
	free(_ax);  free(_r);  free(_tmps1); free(_tmps2); free(_vec_even); free(_vec_odd);
#else
	free(ax); free(r); free(tmps1); free(tmps2); free(vec_even); free(vec_odd);
#endif
	free(evecs); free(evals); free(H); free(HU); free(Hinv);
	free(initwork); free(tmpv1); free(zheev_work);
	free(hevals); free(zheev_rwork); free(IPIV);free(ef_evs);

	return 0;
}


#include <stdlib.h>
#include <complex.h>
//#include <qhg.h>
//#include <qhg_types.h>
//#include <qhg_eigen_field_linalg.h>
#include <mpi.h>
#include <math.h>
#include"solver/eigen_topological.h"
#include <stdio.h>
#include <stdlib.h>

_Complex double top_suscept_all_eigvcs( eigen_field* eigf)
{
    return top_suscept_subset( eigf,eigf->nevs);
}

//change this so subset calls the main one with a new structure
_Complex double top_suscept_subset( eigen_field* eigf,int max_nevs){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int eo_multiple =2; //1 was used for testing purposes when only even-odd e-vectors available

    int saved_nevs =eigf->nevs;
   // eigen_field eigf =ef;
    if(max_nevs<eigf->nevs)
    {
        eigf->nevs= max_nevs;
    }

    _Complex double* local_eig_product =local_eigvecs_g5_self_in_prod(eigf,eo_multiple);

    _Complex double local_trace_sqr=0 +0*_Complex_I;
    _Complex double local_trace=0 + 0*_Complex_I;

    for(int i=0; i<eigf->nevs;i++)
    {
        if(isnan(creal(*(local_eig_product+i))) ==1)
        {
            printf("Issue is eig_product");
        }
  //      if(isnan(creal(*(eigf->evls+i)))==1)
        if(isnan((*(eigf->evls+i)))==1)
        {
            printf("Issue is evls");
        }
        local_trace_sqr += *(local_eig_product+i)/((*(eigf->evls+i))*(*(eigf->evls+i)));
        local_trace += *(local_eig_product+i)/(*(eigf->evls+i));
    }
/* printf("The loc trace sqr is %f + %fi on rank %d\n",creal(local_trace_sqr),cimag(local_trace_sqr),rank);
 *    printf("The loc trace is %f + %fi on rank %d\n",creal(local_trace),cimag(local_trace),rank); */
    _Complex double global_trace_sqr=0;
    _Complex double global_trace=0;
    MPI_Allreduce(&local_trace_sqr, &global_trace_sqr, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_trace,&global_trace, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    eigf->nevs = saved_nevs;
    free(local_eig_product);
    return global_trace_sqr*global_trace;
}


_Complex double top_suscept_product_print(eigen_field*  eigf,int max_nevs,char *v_g5_logfile){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int eo_multiple =2; //1 was used for testing purposes when only even-odd e-vectors available

    //eigen_field eigf = eigf;

    int saved_nevs =eigf->nevs;
    if(max_nevs<eigf->nevs)
    {
        eigf->nevs= max_nevs;
    }

    _Complex double* local_eig_product =local_eigvecs_g5_self_in_prod(eigf,eo_multiple);

    _Complex double* global_eig_product =malloc(sizeof(_Complex double)*eigf->nevs);


    if(global_eig_product==NULL)
    {
	if(rank == 0){
	      fprintf(stderr,"Ran out of memory inside top_suscept_product_print.\n");
        }
	exit(1);
    }

    for(int i=0; i<eigf->nevs;i++)
    {
        *global_eig_product =0;
    }

    MPI_Allreduce(local_eig_product,global_eig_product, 2*(eigf->nevs), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(rank==0)
    {
        FILE * fp;
	printf("The eigenvalue logfile is %s | ", v_g5_logfile);
     //   fp = fopen ("v_g5_product.txt", "w+");
        fp = fopen (v_g5_logfile, "w+");
	

        for(int i=0;i<eigf->nevs;i++)
        {
            fprintf(fp,"%d %.10e %.10e Eigenvalue %.10e \n",i,creal(*(global_eig_product+i)),cimag(*(global_eig_product+i)),*(eigf->evls+i));
//            printf("Product %d %.10e %.10e Eigenvalue %.10e \n",i,creal(*(global_eig_product+i)),cimag(*(global_eig_product+i)),*(eigf->evls+i));
            printf("Product %d %.10e %.10e Eigenvalue %.10e \n",i,creal(*(global_eig_product+i)),cimag(*(global_eig_product+i)),*(eigf->evls+i));
        }

        fclose(fp);
    }
    free(local_eig_product);
    free(global_eig_product);
    eigf->nevs = saved_nevs;
    return 1;
}

void eigen_field_init( _Complex double*  evcs, double* evls,int nevs,int lvol,eigen_field* ef)
{
  ef->evcs = evcs;
  ef->evls =  evls;
 // ef->lat = lat;
  ef->lvol=lvol;
  ef->nevs = nevs;
}

_Complex double* local_eigvecs_g5_self_in_prod(const eigen_field* const eigf,int eo_multiple){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    const int eigvec_size = 12*(eigf->lvol/2)*eo_multiple;
    _Complex double* l_norm_array = malloc(sizeof(_Complex double)*eigf->nevs);
    if(l_norm_array==NULL)
    {
	if(rank == 0){
	      fprintf(stderr,"Ran out of memory inside eigenvector inner product.\n");
        }
	exit(1);
    }
	
    for(int i =0; i<eigf->nevs; i++)
    {
        _Complex double* ptr = (eigf->evcs)+i*eigvec_size;
        l_norm_array[i] = local_g5_self_in_prod2(ptr,eigvec_size);
//       printf("The self value is %f + %fi \n, index %d  on rank %d",creal(*(l_norm_array+i)),cimag(*(l_norm_array+i)),i,rank);
       /*printf("The eigvalue is %f + %fi \n",creal(*(eigf->evls+i)),cimag(*(l_norm_array+i)));*/

    }
//      printf("The self value is %f + %fi",creal(*l_norm_array),cimag(*l_norm_array));
    return l_norm_array;
}


_Complex double local_g5_self_in_prod(const _Complex double* const evcs, const int N)
{
    _Complex double res;
    double real_part=0;
    double imag_part=0;
   #pragma omp parallel reduction(+:real_part,imag_part)
    {
        _Complex double ks,kc,ds,tr,ts,tt;
        const _Complex double *s;

        ks = 0.0;
        kc = 0.0;
//12 is from NS*NC, KAHAN ALGORITHM USED FOR SUMMATION
        #pragma omp for
        for (int ix  =  0; ix < N; ix+=12) {
            s = evcs + ix;
            ds = conj(*s) * (*s) +
                conj(*(s+1)) * *(s+1) +
                conj(*(s+2)) * *(s+2) +
                conj(*(s+3)) * *(s+3) +
                conj(*(s+4)) * *(s+4) +
                conj(*(s+5)) * *(s+5) -
                conj(*(s+6)) * *(s+6) -
                conj(*(s+7)) * *(s+7) -
                conj(*(s+8)) * *(s+8) -
                conj(*(s+9)) * *(s+9) -
                conj(*(s+10)) * *(s+10) -
                conj(*(s+11)) * *(s+11);

            tr = ds + kc;
            ts = tr + ks;
            tt = ts-ks;
            ks = ts;
            kc = tr-tt;
            //printf("This is ds, %.10e + %.10eI",creal(ds),cimag(ds));
        }
        kc=ks+kc;

        real_part += creal(kc);
        imag_part += cimag(kc);
    }
    return real_part + _Complex_I*imag_part;
}





_Complex double* local_eigvecs_g5_in_prod(const eigen_field* const eigf,int eo_multiple){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	const int eigvec_size = 12*(eigf->lvol/2)*eo_multiple;
	_Complex double* C_top_matrix = malloc(sizeof(_Complex double)*eigf->nevs*eigf->nevs);


	if(C_top_matrix==NULL)
	{
		if(rank == 0){
			fprintf(stderr,"Ran out of memory inside eigenvector inner product.\n");
		}
		exit(1);
	}

	for(int i =0; i<eigf->nevs; i++)
	{

		_Complex double* ptr_i = (eigf->evcs)+i*eigvec_size;
		for(int j =0; j<eigf->nevs; j++)
		{
			_Complex double* ptr_j = (eigf->evcs)+j*eigvec_size;
			C_top_matrix[(eigf->nevs)*i +j] = local_g5_in_prod(ptr_i,ptr_j,eigvec_size);
			//			printf("The self value is %f + %fi \n, index %d  on rank %d",creal(*(l_norm_array+i)),cimag(*(l_norm_array+i)),i,rank);
			//   printf("The eigvalue is %f + %fi \n",creal(*(eigf->evls+i)),cimag(*(l_norm_array+i));
		}
	}
	/* printf("The self value is %f + %fi",creal(*l_norm_array),cimag(*l_norm_array));*/
	return C_top_matrix;
}

_Complex double Q_squared_1( eigen_field*  eigf,int max_nevs){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int eo_multiple =2; /*1 was used for testing purposes when only even-odd e-vectors available*/

	int saved_nevs =eigf->nevs;
	//eigen_field eigf=ef;
	if(max_nevs<eigf->nevs)
	{
		eigf->nevs= max_nevs;
	}

	_Complex double* local_eig_product =local_eigvecs_g5_in_prod(eigf,eo_multiple);
	_Complex double* global_i_j_eig_prod = malloc(sizeof(_Complex double)*eigf->nevs*eigf->nevs);

	if(global_i_j_eig_prod==NULL)
	{
		if(rank == 0){
			fprintf(stderr,"Ran out of memory inside eigenvector inner product.\n");
		}
		exit(1);
	}

	MPI_Allreduce(local_eig_product, global_i_j_eig_prod, 2*eigf->nevs*eigf->nevs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	/*  MPI_Allreduce(&local_trace,&global_trace, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);*/

	_Complex double full_sum=0 +0*_Complex_I;
	_Complex double diag_sum=0 +0*_Complex_I;
	for(int i =0;i<eigf->nevs;++i)
	{
		_Complex double* ptr = global_i_j_eig_prod+i*eigf->nevs+ i;
		//full_sum+= conj(*ptr)*(*ptr);
		diag_sum+=*ptr;
	}

	eigf->nevs=saved_nevs;
	free(local_eig_product);
	free(global_i_j_eig_prod);
	return (diag_sum*diag_sum);
}


_Complex double Q_squared_2(eigen_field*  eigf,int max_nevs){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int eo_multiple =2; /*1 was used for testing purposes when only even-odd e-vectors available*/

	int saved_nevs =eigf->nevs;
//	eigen_field eigf=ef;
	
	if(max_nevs<eigf->nevs)
	{
		eigf->nevs= max_nevs;
	}

	_Complex double* local_eig_product =local_eigvecs_g5_self_in_prod(eigf,eo_multiple);
	_Complex double* global_eig_product = malloc(sizeof(_Complex double)*eigf->nevs);

	if(global_eig_product==NULL)
	{
		if(rank == 0){
			fprintf(stderr,"Ran out of memory inside eigenvector inner product.\n");
		}
		exit(1);
	}

	MPI_Allreduce(local_eig_product, global_eig_product, 2*eigf->nevs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	/*  MPI_Allreduce(&local_trace,&global_trace, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);*/

	_Complex double diag_sum=0 +0*_Complex_I;

	for(int i =0;i<eigf->nevs;++i)
	{
		_Complex double* ptr = global_eig_product+i;
		diag_sum+=*ptr;
	}
	
	eigf->nevs =saved_nevs;
	free(local_eig_product);
	free(global_eig_product);
	return (diag_sum*diag_sum);
}
/*
_Complex double renorm_top_suscept_alt(eigen_field* eigf,int max_nevs){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int eo_multiple =2; //1 was used for testing purposes when only even-odd e-vectors available

//	eigen_field eigf=ef;
	int saved_nevs =eigf->nevs;
	if(max_nevs<eigf->nevs)
	{
		eigf->nevs= max_nevs;
	}

	_Complex double* local_eig_product =local_eigvecs_g5_in_prod(eigf,eo_multiple);
	_Complex double* global_i_j_eig_prod = malloc(sizeof(_Complex double)*eigf->nevs*eigf->nevs);

	if(global_i_j_eig_prod==NULL)
	{
		if(rank == 0){
			fprintf(stderr,"Ran out of memory inside eigenvector inner product.\n");
		}
		exit(1);
	}

	MPI_Allreduce(local_eig_product, global_i_j_eig_prod, 2*eigf->nevs*eigf->nevs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// MPI_Allreduce(&local_trace,&global_trace, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	_Complex double full_sum=0 +0*_Complex_I;
	_Complex double diag_sum=0 +0*_Complex_I;
	for(int i =0;i<eigf->nevs;++i)
	{
		for(int j=0;j<eigf->nevs;++j){
			_Complex double* ptr = global_i_j_eig_prod+i*eigf->nevs + j;
			full_sum+= conj(*ptr)*(*ptr);
			if(i==j){
				diag_sum+=*ptr;
			}
		}
	}

	eigf->nevs=saved_nevs;
	free(local_eig_product);
	free(global_i_j_eig_prod);
	return eigf->nevs*(diag_sum*diag_sum)/full_sum;
}

*/

_Complex double local_g5_self_in_prod2(const _Complex double* const evcs, const int N){

	return local_g5_in_prod(evcs, evcs, N);
}


_Complex double local_g5_in_prod(const _Complex double* const evcs_i,const _Complex double* const evcs_j, const int N)
{
    _Complex double res;
    double real_part=0;
    double imag_part=0;

   #pragma omp parallel reduction(+:real_part,imag_part)
    {
        _Complex double ks,kc,ds,tr,ts,tt;
        const _Complex double *s;
        const _Complex double *s2;
        ks = 0.0;
        kc = 0.0;

    #pragma omp for
        for (int ix  =  0; ix < N; ix+=12) {//12 is size of spinor, this is a VERY fixed number in this program.
            s = evcs_i + ix;
            s2 = evcs_j + ix;

            ds = conj(*s) * (*s2) +
                conj(*(s+1)) * *(s2+1) +
                conj(*(s+2)) * *(s2+2) +
                conj(*(s+3)) * *(s2+3) +
                conj(*(s+4)) * *(s2+4) +
                conj(*(s+5)) * *(s2+5) -
                conj(*(s+6)) * *(s2+6) -
                conj(*(s+7)) * *(s2+7) -
                conj(*(s+8)) * *(s2+8) -
                conj(*(s+9)) * *(s2+9) -
                conj(*(s+10)) * *(s2+10) -
                conj(*(s+11)) * *(s2+11);

            tr = ds + kc;
            ts = tr + ks;
            tt = ts-ks;
            ks = ts;
            kc = tr-tt;
        }
        kc=ks+kc;
        real_part += creal(kc);
        imag_part += cimag(kc);
    }
    return real_part + _Complex_I*imag_part;
}

_Complex double v_g5_product_compute(_Complex double* evec,int lvol){

    eigen_field eigf;
    eigen_field_init(evec,NULL,1,lvol,&eigf);   

    _Complex double* local_eig_product =local_eigvecs_g5_self_in_prod(&eigf,2);

    _Complex double global_eig_product =0;

    MPI_Allreduce(local_eig_product,&global_eig_product, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    free(local_eig_product);
    return global_eig_product;

}
void v_g5_product_print(_Complex double* v_g5_product,double* evals, int nconv,char* v_g5_logfile){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
    if(rank==0)
    {
        FILE * fp;
	printf("The eigenvalue logfile is %s | ", v_g5_logfile);
     //   fp = fopen ("v_g5_product.txt", "w+");
        fp = fopen (v_g5_logfile, "w+");
	

        for(int i=0;i<nconv;i++)
        {
            fprintf(fp,"%d %.10e %.10e Eigenvalue %.10e \n",i,creal(*(v_g5_product+i)),cimag(*(v_g5_product+i)),*(evals+i));
            printf("V_G5_Product: %d %.10e %.10e Eigenvalue %.10e \n",i,creal(*(v_g5_product+i)),cimag(*(v_g5_product+i)),*(evals+i));
        }

        fclose(fp);
    }
}



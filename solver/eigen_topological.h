
#ifndef TOPOLOGIC_H
#define TOPOLOGIC_H 

#include <complex.h>
//#include <stdio.h>

typedef struct {
	_Complex double *evcs;
	double *evls;
//	qhg_lattice *lat;
	int lvol;//local_volume
    int nevs;
} eigen_field;

//#include <complex.h>
//#include <qhg_eigen_field.h>
_Complex double top_suscept_subset(  eigen_field*   eigf,int max_nevs);
_Complex double top_suscept_product_print( eigen_field*  eigf,int max_nevs);

void eigen_field_init(_Complex double*  evcs,double* evls,int nevs,int lvol,eigen_field* ef);

_Complex double local_g5_self_in_prod(const _Complex double* const evcs, const int N);
_Complex double local_g5_self_in_prod2(const _Complex double* const evcs, const int N);
_Complex double* local_eigvecs_g5_self_in_prod(const eigen_field* const eigf,int eo_multiple);
_Complex double* local_eigvecs_g5_in_prod(const eigen_field* const eigf,int eo_multiple);
_Complex double local_g5_in_prod(const _Complex double* const evcs_i,const _Complex double* const evcs_j, const int N);
_Complex double Q_squared_1( eigen_field*  eigf,int max_nevs);
_Complex double Q_squared_2( eigen_field* eigf,int max_nevs);
//_Complex double renorm_top_suscept_alt( eigen_field*  eigf,int max_nevs);

#endif

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>

#define C    2.9979e10
#define PI   (2.0*acos(0.0))

#define wave1 (5006.84-4958.91)
#define wave2 (5006.84-4861.34)
#define  flux1  (1.0/3.0)

extern const int nlist_max, ncode_max;
extern char **list;
extern double *o3flux, *o3center;
extern int nlist, ncode;

extern const double o3flux_std;  // in a unit of 1.0e-14 erg/s

extern double **hbb, **o3b, **o3n, **hbn, **optflux, **optslope;
extern double *chi2, *chi2opt;

extern double *date;

extern char **code;
extern int *obs_num, *code_idx;

extern double *ps_scale, *es_scale, *ps_scale_err, *es_scale_err;;

//
extern double * Smat, * Nmat, * INmat, * ISmat, * Qmat, * IQmat, * eigens, * eigens_vecs, * eigens_mat, 
       * Cmat, * ICmat, * N0mat;
extern double * Fcon, * Fcon_err, * Fhb, *Fhb_err;

extern int n_mcmc, nbuilt;
extern int const ntheta_max;
extern double *cov_matrix;

//GSL
/* for GSL */
extern const gsl_rng_type * gsl_T;
extern gsl_rng * gsl_r;

//reconstruct
extern int n_recon;
extern double * USmat;
extern double * date_recon, * Fcon_recon, * Fcon_err_recon, * Fhb_recon, * Fhb_err_recon;
extern double * var_con_best, * var_con_best_err, * var_hb_best, * var_hb_best_err;

#endif

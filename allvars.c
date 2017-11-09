#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"

const int nd_max=2000, ncode_max=100;
char **list;
double *o3flux, *o3center;
int ncode, nd_cont, nd_line;

const double o3flux_std = 1.0; //(5.58e-13/1.0e-14);  // in a unit of 1.0e-14 erg/s

double **hbb, **o3b, **o3n, **hbn, **optflux, **optslope;
double *chi2, *chi2opt;
double *date_cont, *date_line;

char **code;
int *obs_num_cont, *obs_num_line, *code_idx_cont, *code_idx_line;
double *ps_scale, *es_scale, *ps_scale_err, *es_scale_err;


// MCMC
double * Smat, * Nmat, * INmat, * ISmat, * Qmat, * IQmat, * eigens, * eigens_vecs, * eigens_mat, 
       * Cmat, * ICmat, * N0mat;
double * Fcon, * Fcon_err, * Fhb, *Fhb_err;

int n_mcmc, nbuilt;
int const ntheta_max=50;
double *cov_matrix;

//GSL
/* for GSL */
const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;


//reconstruct
int n_recon;
double * USmat;
double * date_recon, * Fcon_recon, * Fcon_err_recon, * Fhb_recon, * Fhb_err_recon;
double * var_hb_best, * var_hb_best_err, * var_con_best, * var_con_best_err;



PARSET parset;
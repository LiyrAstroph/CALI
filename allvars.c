#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"

const int nd_max=2000, ncode_max=100;

int ncode, nd_cont, nd_line, ndrw;

double **hbb, **optflux;

double date_span_cont, date_span_line;

double *date_cont, *date_line;
double *date_cont_org, *date_line_org, **optflux_org, **hbb_org;
size_t *perm_cont, *perm_line;


char **code;
int *obs_num_cont, *obs_num_line, *code_idx_cont, *code_idx_line;
int *code_idx_cont_org, *code_idx_line_org;

double *ps_scale, *es_scale, *ps_scale_err, *es_scale_err;


// MCMC
double * Smat, * Nmat, * INmat, * ISmat, * Qmat, * IQmat, 
       * Cmat, * ICmat, * N0mat;
double * Fcon, * Fcon_err, * Fhb, *Fhb_err;

double *workspace;

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
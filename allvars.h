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

extern const int nd_max, ncode_max;
extern char **list;
extern double *o3flux, *o3center;
extern int ncode, nd_cont, nd_line;

extern const double o3flux_std;  // in a unit of 1.0e-14 erg/s

extern double **hbb, **o3b, **o3n, **hbn, **optflux, **optslope;
extern double *chi2, *chi2opt;

extern double *date_cont, *date_line;

extern char **code;
extern int *obs_num_cont, *obs_num_line, *code_idx_cont, *code_idx_line;

extern double *ps_scale, *es_scale, *ps_scale_err, *es_scale_err;;

//
extern double * Smat, * Nmat, * INmat, * ISmat, * Qmat, * IQmat, 
       * Cmat, * ICmat, * N0mat;
extern double * Fcon, * Fcon_err, * Fhb, *Fhb_err;

extern double *workspace;

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


typedef struct
{
  char file_param[256];
  char file_cont[256];
  char file_line[256];

  int n_mcmc, n_builtin;

  int flag_line;

}PARSET;
extern PARSET parset;

#endif

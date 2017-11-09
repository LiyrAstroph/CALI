#include <gsl/gsl_histogram.h>


//prototype

int output_optical(char *fname);
int output_hb(char *fname);

int memory_malloc();
int calibrate();

//mcmc
int mcmc();
int mcmc_pt();
void mcmc_stats();
void mcmc_memory_init();
double prob_variability(double *var_con, double *var_hb);
void set_covar_mat_con(double sigma, double tau, double alpha);
void set_covar_mat_hb(double sigma, double tau, double alpha);
void scale_flux();
void set_scale_err(double *theta_var);
void set_scale(double *theta);
void set_mcmc_param(double *theta, double *stepsize, double (*range)[2], int ntheta);
void init_cov_matrix(double *stepsize, int ntheta);
void get_cov_matrix(double *theta, int nstep, int ntheta);
void get_cov_matrix_diag(double *theta, int nstep, int ntheta);
int read_dataset();
double prob_variability_beta(double *var_con, double *var_hb, double beta);

/* reconstruct */
void reconstruct_int();
int reconstruct_con(double *var_con);
int reconstruct_hb(double *var_hb);
void set_covar_Umat_con(double sigma, double tau, double alpha);
void set_covar_Umat_hb(double sigma, double tau, double alpha);


void read_parset();

/* matrix operations */
void inverse_mat(double *a, int n, int *info);
double det_mat(double *a, int n, int *info);
double lndet_mat(double *a, int n, int *info);
void display_mat(double *a, int m, int n);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_tranposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);
int get_histogram_val_error(gsl_histogram *hist, double *val, double *err_lo, double *err_hi);

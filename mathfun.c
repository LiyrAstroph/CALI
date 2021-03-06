#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cblas.h>
//#include <clapack.h>
#include <lapacke.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "proto.h"

void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, k);
}
void multiply_mat_MN_tranposeA(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, k);
}

void multiply_mat(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}
void multiply_mat_transposeA(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

void multiply_mat_transposeB(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

void multiply_matvec_transposeA(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

void multiply_matvec(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/* y(m) = a(m, n) * x(n) */
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/* C(m*n) = A^T(m*k) * B(k*n) */
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1.0f
                             , a, m, b, n, 0.0f, c, n);
}
/* C(m*n) = A(m*k) * B^T(k*n) */
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0f
                             , a, k, b, k, 0.0f, c, n);
}

void inverse_mat(double * a, int n, int *info)
{
  int * ipiv;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);

  free(ipiv);
  return;
}

void eigen_sym_mat(double *a, int n, double *val, int *info)
{
    char jobz='V', uplo='U';

/* store the eigenvectors  in a by rows.
 * store the eigenvalues in val in ascending order.
 */
//    dsyev_(&jobz, &uplo, &n, a, &n, val, work, &lwork, info);

/* store the eigenvectors  in a by columns.
 * store the eigenvalues in val in ascending order.
 */
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, val);
    return;
}

void multiply_vec2mat(double * x, double * a, int n)
{
//  cblas_dsyr(CblasRowMajor, CblasUpper, n, 1.0f, x, 1, a, n);
  int i, j;
  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      a[i*n+j] = a[j*n+i] = x[i]*x[j];
    }
}

/* get determinant of matrix A */
double det_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double det;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("Wrong!\n");
    exit(-1);
  }

  det = 1.0;
  for(i=0; i<n; i++)
  {
    det *= a[i*n+i];
    if (ipiv[i] != i)
    {
      ipiv[ipiv[i]] = ipiv[i];
      det = -det;
    }
  }
  det=fabs(det);
  free(ipiv);
/*  char uplo='U';
  int i;
  double det;
//  dpotrf_(&uplo, &n, a, &n, info);
  LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  det=1.0;
  for(i=0;i<n;i++)
    det *= a[i*n+i];
  det *=det;*/

/*  int s;
  double det;
  gsl_matrix_view m = gsl_matrix_view_array (a, n, n);
  gsl_permutation * p = gsl_permutation_alloc (n);
  gsl_linalg_LU_decomp(&m.matrix, p, &s );
  det = gsl_linalg_LU_det(&m.matrix, s);*/

  return det;
}

double lndet_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double lndet;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("Wrong!\n");
    exit(-1);
  }

  lndet = 0.0;
  for(i=0; i<n; i++)
  {
    lndet += log(a[i*n+i]);
  }
  free(ipiv);
  return lndet;
}

void Chol_decomp_U(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'U';
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
    return;
  }
  else if (*info>0)
  {
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
    return;
  }
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      a[i*n+j] = 0.0;
  return;
}

/*
 * semiseparable matrix
 */
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr, double *W, double *D, double *phi)
{
  int i;
  double S, A;
  phi[0] = 0.0;
  for(i=1; i<n; i++)
  {
    phi[i] = exp(-c1 * (t[i] - t[i-1]));
  }

  S = 0.0;
  A = sigma[0]*sigma[0] + syserr*syserr + a1;
  D[0] = A;
  W[0] = 1.0/D[0];
  for(i=1; i<n; i++)
  {
    S = phi[i]*phi[i] * (S + D[i-1]*W[i-1]*W[i-1]);
    A = sigma[i]*sigma[i] + syserr*syserr + a1;
    D[i] = A - a1 * a1 * S;
    W[i] = 1.0/D[i] * (1.0 - a1*S);
  }
}
/*
 * z = C^-1 x y
 *
 * y is a vector
 */
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z)
{
  int i;
  double f, g;

  // forward substitution
  f = 0.0;
  z[0] = y[0];
  for(i=1; i<n;i++)
  {
    f = phi[i] * (f + W[i-1] * z[i-1]);
    z[i] = y[i] - a1*f;
  }

  //backward substitution
  g = 0.0;
  z[n-1] = z[n-1]/D[n-1];
  for(i=n-2; i>=0; i--)
  {
    g = phi[i+1] *(g + a1*z[i+1]);
    z[i] = z[i]/D[i] - W[i]*g;
  }
}
/*
 * Z = C^-1 x Y
 * 
 * Y is an (nxm) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0*m+j];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i*m+j] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

/*
 * Z = C^-1 x Y^T
 * 
 * Y is an (mxn) matrix. 
 * Note that Y is row-major
 */
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z)
{
  int i, j;
  double f, g;

  // forward substitution
  for(j=0; j<m; j++)
  {
    f = 0.0;
    Z[0*m+j] = Y[0+j*n];
    for(i=1; i<n;i++)
    {
      f = phi[i] * (f + W[i-1] * Z[(i-1)*m + j]);
      Z[i*m+j] = Y[i+j*n] - a1*f;
    }
  }

  //backward substitution
  for(j=0; j<m; j++)
  {
    g = 0.0;
    Z[(n-1)*m+j] = Z[(n-1)*m+j]/D[n-1];
    for(i=n-2; i>=0; i--)
    {
      g = phi[i+1] *(g + a1*Z[(i+1)*m+j]);
      Z[i*m+j] = Z[i*m+j]/D[i] - W[i]*g;
    }
  }
}

void display_mat(double *a, int m, int n)
{
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%e\t", a[i*n+j]);
        }
        printf("\n");
    }
}

double ** matrix_malloc(int n1, int n2)
{
  double ** mat;
  int i;

  if(!(mat = malloc(n1*sizeof(double*))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  for(i=0; i<n1; i++)
  {
    if(!(mat[i] = malloc(n2*sizeof(double))))
    {
      fprintf(stderr, "Unable to allocate the matrix!\n");
      exit(-1);
    }
  }
  return mat;
}

double * array_malloc(int n)
{
  double *array;

  if(!(array = malloc(n*sizeof(double))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  return array;
}

/*
 * get the best estimate of parameter and its errors from GLS histogram
 */
int get_histogram_val_error(gsl_histogram *hist, double *val, double *err_lo, double *err_hi)
{
  int idx, i, nh;
  double max;

  nh = gsl_histogram_bins(hist);
  idx = gsl_histogram_max_bin(hist);
  max = gsl_histogram_max_val(hist);
  *val = 0.5*(hist->range[idx] + hist->range[idx+1]);

  max *= exp(-0.5);

  for(i=idx-1; i>0; i--)
  {
    if(hist->bin[i]- max >= 0 && hist->bin[i-1] - max <=0 )
    {
      *err_lo = 0.25 * ( hist->range[i-1] + hist->range[i+1] + 2.0*hist->range[i] );
      break;
    }
  } 

  if(i==0)
    *err_lo = hist->range[0];

  for(i=idx; i<nh-1; i++)
  {
    if(hist->bin[i]-max >= 0 && hist->bin[i+1] - max <=0 )
    {
      *err_hi = 0.25 * ( hist->range[i] + hist->range[i+2] + 2.0*hist->range[i+1] );
      break;
    }
  }  
  if(i==nh-1)
    *err_hi = hist->range[nh];
  
  return 0;
}

void get_cov_matrix(double *theta, int nstep, int ntheta)
{
  int i, j;
  double cov;
  
  for(i=0; i<ntheta; i++)
  {
    cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[i*nstep], 1, nstep);
    cov = cov>1.0e-8? cov : 1.0e-8;
    cov_matrix[i*ntheta + i] = cov;
    for(j=0; j<i; j++)
    {
      cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[j*nstep], 1, nstep);
      cov = cov>1.0e-8?cov:0.0;
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = cov;
    }
  }
}

void get_cov_matrix_diag(double *theta, int nstep, int ntheta)
{
  int i, j;
  double cov;
  
  for(i=0; i<ntheta; i++)
  {
    cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[i*nstep], 1, nstep);
    cov = cov>1.0e-6 ? cov : 1.0e-6;
    if(i<4)cov = cov<0.1 ? cov : 0.1;
    if(i>=4)cov = cov<1.0e-4 ? cov : 1.0e-4;
    cov_matrix[i*ntheta + i] = cov;
    for(j=0; j<i; j++)
    {
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = 0.0;
    }
  } 
/*  for(i=0; i<4; i++)
    for(j=0; j<i; j++)
    {
      cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[j*nstep], 1, nstep);
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = cov;
    }*/
}

double mod(double y, double x)
{
  if(x > 0.0)
  {
    return (y/x - floor(y/x))*x;
  }
  else if(x == 0.0)
  {
    return 0.0;
  }
  else
  {
    printf("Warning in mod(double, double) %e\n", x);
    exit(0);
  }
  
}

void wrap(double *x, double min, double max)
{
  *x = mod(*x - min, max - min) + min;
}

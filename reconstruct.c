#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "proto.h"

int reconstruct_con(double *var_con)
{
  FILE *fout;

  double sigma, tau, alpha, ave_con, lambda, errtemp;
  double *Larr, *ybuf, *y, *Lmat, *Vmat, *Sstar;
  int i, info;

  y = array_malloc(nlist);
  ybuf = array_malloc(nlist);
  Larr = array_malloc(nlist);
  Sstar = array_malloc(nlist);
  Lmat = array_malloc(nlist*nlist);
  Vmat = array_malloc(nlist*nlist);

  sigma = var_con[0];
  tau = var_con[1];
  alpha = var_con[2];
  
  printf("%f\t%f\t%f\n", sigma, tau, alpha);

  set_covar_mat_con(sigma, tau, alpha);
  set_covar_Umat(sigma, tau, alpha);

  for(i=0;i<nlist*nlist; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }

  memcpy(ICmat, Cmat, nlist*nlist*sizeof(double));
  inverse_mat(ICmat, nlist, &info);

  for(i=0;i<nlist;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nlist, ybuf);
  lambda = 0.0;
  for(i=0;i<nlist;i++)lambda+=Larr[i]*ybuf[i];
  multiply_matvec(ICmat, Fcon, nlist, ybuf);
  ave_con = 0.0;
  for(i=0;i<nlist;i++)ave_con+=Larr[i]*ybuf[i];
  ave_con /=lambda;

  for(i=0; i<nlist; i++)y[i] = Fcon[i] - ave_con;
  multiply_matvec(ICmat, y, nlist, ybuf);
  multiply_matvec_MN(USmat, n_recon, nlist, ybuf, Fcon_recon);
  for(i=0; i<n_recon; i++)Fcon_recon[i] = Fcon_recon[i] + ave_con;

  /* variance estimate */
  for(i=0;i<nlist*nlist; i++)Lmat[i]=1.0;
  multiply_vec2mat(Larr, Lmat, nlist);

  multiply_mat(Lmat, ICmat, Vmat, nlist);
  multiply_mat(ICmat, Vmat, Lmat, nlist);
//!  display_mat(Lmat, n_con_data, n_con_data);
  for(i=0; i<nlist*nlist; i++)Vmat[i] = ICmat[i];// - Lmat[i]/lambda;
  for(i=0; i<n_recon; i++)
  {
    memcpy(Sstar, &USmat[i*nlist], nlist*sizeof(double));
    multiply_matvec(Vmat, Sstar, nlist, ybuf);
    Fcon_err_recon[i] = sigma*sigma - cblas_ddot(nlist, Sstar, 1, ybuf, 1);
    errtemp=cblas_ddot(nlist, Larr, 1, ybuf, 1)-1.0;
    Fcon_err_recon[i] += errtemp * errtemp / lambda;
//    printf("%f\n", errtemp*errtemp/lambda);
    Fcon_err_recon[i] = sqrt(Fcon_err_recon[i]);
  }

  fout = fopen("opt_recon.txt", "w");

  for(i=0; i<n_recon; i++)
    fprintf(fout, "%f\t%f\t%f\n", date_recon[i], Fcon_recon[i], Fcon_err_recon[i]);

  fclose(fout);
  return 0;
}

int reconstruct_hb(double *var_hb)
{
  FILE *fout;

  double sigma, tau, alpha, ave_con, lambda, errtemp;
  double *Larr, *ybuf, *y, *Lmat, *Vmat, *Sstar;
  int i, info;

  y = array_malloc(nlist);
  ybuf = array_malloc(nlist);
  Larr = array_malloc(nlist);
  Sstar = array_malloc(nlist);
  Lmat = array_malloc(nlist*nlist);
  Vmat = array_malloc(nlist*nlist);

  sigma = var_hb[0];
  tau = var_hb[1];
  alpha = var_hb[2];
  
  printf("%f\t%f\t%f\n", sigma, tau, alpha);

  set_covar_mat_hb(sigma, tau, alpha);
  set_covar_Umat(sigma, tau, alpha);

  for(i=0;i<nlist*nlist; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }

  memcpy(ICmat, Cmat, nlist*nlist*sizeof(double));
  inverse_mat(ICmat, nlist, &info);

  for(i=0;i<nlist;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nlist, ybuf);
  lambda = 0.0;
  for(i=0;i<nlist;i++)lambda+=Larr[i]*ybuf[i];
  multiply_matvec(ICmat, Fhb, nlist, ybuf);
  ave_con = 0.0;
  for(i=0;i<nlist;i++)ave_con+=Larr[i]*ybuf[i];
  ave_con /=lambda;

  for(i=0; i<nlist; i++)y[i] = Fhb[i] - ave_con;
  multiply_matvec(ICmat, y, nlist, ybuf);
  multiply_matvec_MN(USmat, n_recon, nlist, ybuf, Fhb_recon);
  for(i=0; i<n_recon; i++)Fhb_recon[i] = Fhb_recon[i] + ave_con;

  /* variance estimate */
  for(i=0;i<nlist*nlist; i++)Lmat[i]=1.0;
  multiply_vec2mat(Larr, Lmat, nlist);

  multiply_mat(Lmat, ICmat, Vmat, nlist);
  multiply_mat(ICmat, Vmat, Lmat, nlist);
//!  display_mat(Lmat, n_con_data, n_con_data);
  for(i=0; i<nlist*nlist; i++)Vmat[i] = ICmat[i];// - Lmat[i]/lambda;
  for(i=0; i<n_recon; i++)
  {
    memcpy(Sstar, &USmat[i*nlist], nlist*sizeof(double));
    multiply_matvec(Vmat, Sstar, nlist, ybuf);
    Fhb_err_recon[i] = sigma*sigma - cblas_ddot(nlist, Sstar, 1, ybuf, 1);
    errtemp=cblas_ddot(nlist, Larr, 1, ybuf, 1)-1.0;
    Fhb_err_recon[i] += errtemp * errtemp / lambda;
//    printf("%f\n", errtemp*errtemp/lambda);
    Fhb_err_recon[i] = sqrt(Fhb_err_recon[i]);
  }

  fout = fopen("hb_recon.txt", "w");

  for(i=0; i<n_recon; i++)
    fprintf(fout, "%f\t%f\t%f\n", date_recon[i], Fhb_recon[i], Fhb_err_recon[i]);

  fclose(fout);
  return 0;
}

void reconstruct_int()
{
  double date_max, date_min, range, dT;
  int i;

  date_recon = array_malloc(n_recon);
  Fcon_recon = array_malloc(n_recon);
  Fcon_err_recon = array_malloc(n_recon);
  Fhb_recon = array_malloc(n_recon);
  Fhb_err_recon = array_malloc(n_recon);

  USmat = array_malloc(n_recon*nlist);

  date_max = 0.0;
  date_min = 1.0e10;
  for(i=0; i<nlist; i++)
  {
    date_max = fmax(date_max, date[i]);
    date_min = fmin(date_min, date[i]);
  }
  range = date_max - date_min;
  date_max += 0.1*range;
  date_min -= 0.1*range;
  printf("%f\t%f\t%f\n", date_max, date_min, range);

  dT = (date_max - date_min)/(n_recon-1);
  for(i=0; i<n_recon; i++)
  {
    date_recon[i] = date_min + dT * i;
  }

  return;
}

void set_covar_Umat(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n_recon; i++)
  {
    t1 = date_recon[i];
    for(j=0; j<nlist; j++)
    {
      t2 = date[j];
      USmat[i*nlist+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}
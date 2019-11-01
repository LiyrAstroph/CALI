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

  y = array_malloc(nd_cont);
  ybuf = array_malloc(nd_cont);
  Larr = array_malloc(nd_cont);
  Sstar = array_malloc(nd_cont);
  Lmat = array_malloc(nd_cont*nd_cont);
  Vmat = array_malloc(nd_cont*nd_cont);

  tau = var_con[1];
  sigma = var_con[0] * sqrt(tau);
  alpha = var_con[2];
  
  //printf("%f\t%f\t%f\n", sigma, tau, alpha);

  set_covar_mat_con(sigma, tau, alpha);
  set_covar_Umat_con(sigma, tau, alpha);

  for(i=0;i<nd_cont*nd_cont; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }

  memcpy(ICmat, Cmat, nd_cont*nd_cont*sizeof(double));
  inverse_mat(ICmat, nd_cont, &info);

  for(i=0;i<nd_cont;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nd_cont, ybuf);
  lambda = 0.0;
  for(i=0;i<nd_cont;i++)lambda+=Larr[i]*ybuf[i];
  multiply_matvec(ICmat, Fcon, nd_cont, ybuf);
  ave_con = 0.0;
  for(i=0;i<nd_cont;i++)ave_con+=Larr[i]*ybuf[i];
  ave_con /=lambda;

  for(i=0; i<nd_cont; i++)y[i] = Fcon[i] - ave_con;
  multiply_matvec(ICmat, y, nd_cont, ybuf);
  multiply_matvec_MN(USmat, n_recon, nd_cont, ybuf, Fcon_recon);
  for(i=0; i<n_recon; i++)Fcon_recon[i] = Fcon_recon[i] + ave_con;

  /* variance estimate */
  for(i=0;i<nd_cont*nd_cont; i++)Lmat[i]=1.0;
  multiply_vec2mat(Larr, Lmat, nd_cont);

  multiply_mat(Lmat, ICmat, Vmat, nd_cont);
  multiply_mat(ICmat, Vmat, Lmat, nd_cont);
//!  display_mat(Lmat, n_con_data, n_con_data);
  for(i=0; i<nd_cont*nd_cont; i++)Vmat[i] = ICmat[i];// - Lmat[i]/lambda;
  for(i=0; i<n_recon; i++)
  {
    memcpy(Sstar, &USmat[i*nd_cont], nd_cont*sizeof(double));
    multiply_matvec(Vmat, Sstar, nd_cont, ybuf);
    Fcon_err_recon[i] = sigma*sigma - cblas_ddot(nd_cont, Sstar, 1, ybuf, 1);
    errtemp=cblas_ddot(nd_cont, Larr, 1, ybuf, 1)-1.0;
    Fcon_err_recon[i] += errtemp * errtemp / lambda;
//    printf("%f\n", errtemp*errtemp/lambda);
    Fcon_err_recon[i] = sqrt(Fcon_err_recon[i]);
  }

  fout = fopen("cont_recon.txt", "w");

  for(i=0; i<n_recon; i++)
    fprintf(fout, "%f\t%f\t%f\n", date_recon[i], Fcon_recon[i] * flux_mean_cont, Fcon_err_recon[i]* flux_mean_cont);

  fclose(fout);
  return 0;
}

int reconstruct_hb(double *var_hb)
{
  FILE *fout;

  double sigma, tau, alpha, ave_con, lambda, errtemp;
  double *Larr, *ybuf, *y, *Lmat, *Vmat, *Sstar;
  int i, info;

  y = array_malloc(nd_line);
  ybuf = array_malloc(nd_line);
  Larr = array_malloc(nd_line);
  Sstar = array_malloc(nd_line);
  Lmat = array_malloc(nd_line*nd_line);
  Vmat = array_malloc(nd_line*nd_line);

  tau = var_hb[1];
  sigma = var_hb[0]*sqrt(tau);
  alpha = var_hb[2];
  
  //printf("%f\t%f\t%f\n", sigma, tau, alpha);

  set_covar_mat_hb(sigma, tau, alpha);
  set_covar_Umat_hb(sigma, tau, alpha);

  for(i=0;i<nd_line*nd_line; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }

  memcpy(ICmat, Cmat, nd_line*nd_line*sizeof(double));
  inverse_mat(ICmat, nd_line, &info);

  for(i=0;i<nd_line;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nd_line, ybuf);
  lambda = 0.0;
  for(i=0;i<nd_line;i++)lambda+=Larr[i]*ybuf[i];
  multiply_matvec(ICmat, Fhb, nd_line, ybuf);
  ave_con = 0.0;
  for(i=0;i<nd_line;i++)ave_con+=Larr[i]*ybuf[i];
  ave_con /=lambda;

  for(i=0; i<nd_line; i++)y[i] = Fhb[i] - ave_con;
  multiply_matvec(ICmat, y, nd_line, ybuf);
  multiply_matvec_MN(USmat, n_recon, nd_line, ybuf, Fhb_recon);
  for(i=0; i<n_recon; i++)Fhb_recon[i] = Fhb_recon[i] + ave_con;

  /* variance estimate */
  for(i=0;i<nd_line*nd_line; i++)Lmat[i]=1.0;
  multiply_vec2mat(Larr, Lmat, nd_line);

  multiply_mat(Lmat, ICmat, Vmat, nd_line);
  multiply_mat(ICmat, Vmat, Lmat, nd_line);
//!  display_mat(Lmat, n_con_data, n_con_data);
  for(i=0; i<nd_line*nd_line; i++)Vmat[i] = ICmat[i];// - Lmat[i]/lambda;
  for(i=0; i<n_recon; i++)
  {
    memcpy(Sstar, &USmat[i*nd_line], nd_line*sizeof(double));
    multiply_matvec(Vmat, Sstar, nd_line, ybuf);
    Fhb_err_recon[i] = sigma*sigma - cblas_ddot(nd_line, Sstar, 1, ybuf, 1);
    errtemp=cblas_ddot(nd_line, Larr, 1, ybuf, 1)-1.0;
    Fhb_err_recon[i] += errtemp * errtemp / lambda;
//    printf("%f\n", errtemp*errtemp/lambda);
    Fhb_err_recon[i] = sqrt(Fhb_err_recon[i]);
  }

  fout = fopen("line_recon.txt", "w");

  for(i=0; i<n_recon; i++)
    fprintf(fout, "%f\t%f\t%f\n", date_recon[i], Fhb_recon[i] * flux_mean_line, Fhb_err_recon[i]* flux_mean_line);

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
  if(parset.flag_line == 1)
  {
    Fhb_recon = array_malloc(n_recon);
    Fhb_err_recon = array_malloc(n_recon);
  }

  USmat = array_malloc(n_recon*nd_max);

  date_max = 0.0;
  date_min = 1.0e10;
  for(i=0; i<nd_cont; i++)
  {
    date_max = fmax(date_max, date_cont[i]);
    date_min = fmin(date_min, date_cont[i]);
  }

  if(parset.flag_line == 1)
  {
    for(i=0; i<nd_line; i++)
    {
      date_max = fmax(date_max, date_line[i]);
      date_min = fmin(date_min, date_line[i]);
    }
  }

  range = date_max - date_min;
  date_max += 0.1*range;
  date_min -= 0.1*range;
  //printf("%f\t%f\t%f\n", date_max, date_min, range);

  dT = (date_max - date_min)/(n_recon-1);
  for(i=0; i<n_recon; i++)
  {
    date_recon[i] = date_min + dT * i;
  }
  
  return;
}

void reconstruct_end()
{
  free(date_recon);
  free(Fcon_recon);
  free(Fcon_err_recon);

  if(parset.flag_line == 1)
  {
    free(Fhb_recon);
    free(Fhb_err_recon);
  }

  free(USmat);
  return;
}

void set_covar_Umat_con(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n_recon; i++)
  {
    t1 = date_recon[i];
    for(j=0; j<nd_cont; j++)
    {
      t2 = date_cont[j];
      USmat[i*nd_cont+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}

void set_covar_Umat_hb(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n_recon; i++)
  {
    t1 = date_recon[i];
    for(j=0; j<nd_line; j++)
    {
      t2 = date_line[j];
      USmat[i*nd_line+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}
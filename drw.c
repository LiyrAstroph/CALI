#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <cblas.h>
//#include <clapack.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "proto.h"

int mcmc()
{
  FILE *fmcmc_out;

  int const ntheta=ndrw*2+2*(ncode-1), n_cov_update=20000;
  double theta[ntheta], theta_new[ntheta], stepsize[ntheta], range[ntheta][2];  
  double Pmatrix[ntheta*ntheta], Prandvec[ntheta], Py[ntheta];
  double theta_mcmc[ntheta*n_cov_update];
  double scale, var_con[3], var_hb[3];
  double prob, prob_new, sigma, tau, alpha, u, ratio;
  int i, info, flag, istep, iaccpt;
  char fname[100];

  set_mcmc_param(theta, stepsize, range, ntheta);
  
  sprintf(fname, "mcmc.txt");
  fmcmc_out = fopen(fname, "w");
  if(fmcmc_out==NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", fname);
    exit(-1);
  }

  
  init_cov_matrix(stepsize, ntheta);
  memcpy(Pmatrix, cov_matrix, ntheta*ntheta*sizeof(double));
  Chol_decomp_U(Pmatrix, ntheta, &info);

  set_scale(theta);
  scale_flux();
  var_con[0] = exp(theta[0]);
  var_con[1] = exp(theta[1]);
  var_con[3] = 1.0;
  if(parset.flag_line == 1)
  {
    var_hb[0] = exp(theta[2]);
    var_hb[1] = exp(theta[3]);
    var_hb[3] = 1.0;
  }
  prob = prob_variability_semi(var_con, var_hb);
  
  istep =0; 
  iaccpt = 0;
  scale = 0.1;
  
  while (istep < n_mcmc)
  {
   
    for(i=0; i<ntheta; i++)
    {
      Prandvec[i] = gsl_ran_gaussian(gsl_r, 1.0);
    }
    multiply_matvec(Pmatrix, Prandvec, ntheta, Py);
    for(i=0; i<ntheta; i++)
    {
      theta_new[i] = theta[i] + scale  * Py[i];
      wrap(&theta_new[i], range[i][0], range[i][1]);
    }
    
    set_scale(theta_new);
    scale_flux();
    var_con[0] = exp(theta_new[0]);
    var_con[1] = exp(theta_new[1]);
    if(parset.flag_line == 1)
    {
      var_hb[0] = exp(theta_new[2]);
      var_hb[1] = exp(theta_new[3]);
    }
    prob_new = prob_variability_semi(var_con, var_hb);
//    printf("%e %e %e %e\n", prob, prob_new, sigma, tau);
    
    ratio = prob_new - prob;
    if(ratio>0.0)
    {
      memcpy(theta, theta_new, sizeof(theta));
      prob = prob_new;
      iaccpt++;
    }
    else
    {
      u = gsl_rng_uniform(gsl_r);
      if(log(u)<ratio)
      {
        memcpy(theta, theta_new, sizeof(theta));
        prob = prob_new;
        iaccpt++;
      }
    }

    for(i=0; i<ntheta; i++)
        theta_mcmc[i*n_cov_update + istep%n_cov_update] = theta[i];

    if(istep%100==0)printf("%d\n", istep);
    
    fprintf(fmcmc_out, "%d", istep);
    for(i=0;i<ntheta;i++)
    {
      if(i<ndrw*2+ncode-1)
        fprintf(fmcmc_out,"\t%f", theta[i]);
      else
        fprintf(fmcmc_out,"\t%f", theta[i]*flux_mean_cont); // account for flux_mean_cont
    }
    fprintf(fmcmc_out, "\n");

    istep++;

    if(istep%n_cov_update == 0)
    {
      get_cov_matrix_diag(theta_mcmc, n_cov_update, ntheta);
//      display_mat(cov_matrix, ntheta, ntheta);
      //for(i=0; i<ntheta; i++)printf("%f\n", cov_matrix[i*ntheta+i]);
      memcpy(Pmatrix, cov_matrix, ntheta*ntheta*sizeof(double));
      Chol_decomp_U(Pmatrix, ntheta, &info);
    }
  }
  
}

int mcmc_pt()
{

  FILE *fmcmc_out;
  
  int const ntheta=ndrw*2+2*(ncode-1), n_cov_update=20000, nbeta=10;
  
  double const fswap = 1.0/200.0;
  double beta[nbeta], beta_min=1.0e-1, beta_max=1.0, dbeta;
  double ubeta;
  int ibeta;

  double theta[nbeta][ntheta], theta_new[nbeta][ntheta], stepsize[nbeta][ntheta], range[ntheta][2];  
  double Pmatrix[nbeta][ntheta*ntheta], Prandvec[nbeta][ntheta], Py[ntheta];
  double **theta_mcmc, swap[ntheta];

  double scale[nbeta], var_con[3], var_hb[3], var_con1[3], var_hb1[3], var_con2[3], var_hb2[3];
  double prob[nbeta], prob_new[nbeta], sigma, tau, alpha, u, ratio, prob1, prob2;
  int i, j, info, flag, istep, iaccpt[nbeta];
  char fname[100];
  FILE fp[nbeta];

  theta_mcmc = matrix_malloc(nbeta, ntheta*n_cov_update);

  dbeta = (beta_max - beta_min)/(nbeta-1);
  for(j=0; j<nbeta; j++)
  {
    beta[j] = beta_min + j*dbeta;
    //printf("%f\n", beta[j]);
  }

  set_mcmc_param(theta[0], stepsize[0], range, ntheta);
  for(j=1; j<nbeta; j++)
  {
    memcpy(stepsize[j], stepsize[0], sizeof(stepsize[0]));
    memcpy(theta[j], theta[0], sizeof(theta[0]));
  }
  
  sprintf(fname, "mcmc.txt");
  fmcmc_out = fopen(fname, "w");
  if(fmcmc_out==NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", fname);
    exit(-1);
  }

  for(j=0; j<nbeta; j++)
  {
    init_cov_matrix(stepsize[j], ntheta);
    memcpy(Pmatrix[j], cov_matrix, ntheta*ntheta*sizeof(double));
    Chol_decomp_U(Pmatrix[j], ntheta, &info);
  }

  for(j=0; j<nbeta; j++)
  {
    set_scale(theta[j]);
    scale_flux();
    var_con[0] = exp(theta[j][0]);
    var_con[1] = exp(theta[j][1]);
    var_con[2] = 1.0;
    if(parset.flag_line == 1)
    {
      var_hb[0] = exp(theta[j][2]);
      var_hb[1] = exp(theta[j][3]);
      var_hb[2] = 1.0;
    }
    prob[j] = prob_variability_semi_beta(var_con, var_hb, beta[j]);
    scale[j] = 0.1;
    iaccpt[j] = 0;
  }
  
  istep =0; 
  while (istep < n_mcmc)
  {
    for(j=0; j<nbeta; j++)
    {
    
      for(i=0; i<ntheta; i++)
      {
        Prandvec[j][i] = gsl_ran_gaussian(gsl_r, 1.0);
      }
      multiply_matvec(Pmatrix[j], Prandvec[j], ntheta, Py);
      for(i=0; i<ntheta; i++)
      {
        theta_new[j][i] = theta[j][i] + scale[j]   * Py[i];
        wrap(&theta_new[j][i], range[i][0], range[i][1]);
      }

      set_scale(theta_new[j]);
      scale_flux();
      var_con[0] = exp(theta_new[j][0]);
      var_con[1] = exp(theta_new[j][1]);
      if(parset.flag_line == 1)
      {
        var_hb[0] = exp(theta_new[j][2]);
        var_hb[1] = exp(theta_new[j][3]);
      }
      prob_new[j] = prob_variability_semi_beta(var_con, var_hb, beta[j]);

      ratio = prob_new[j]- prob[j] ;
      if(ratio>0.0)
      {
        memcpy(theta[j], theta_new[j], sizeof(theta[j]));
        prob[j] = prob_new[j];
        iaccpt[j]++;
      }
      else
      {
        u = gsl_rng_uniform(gsl_r);
        if(log(u)<ratio)
        {
          memcpy(theta[j], theta_new[j], sizeof(theta[j]));
          prob[j] = prob_new[j];
          iaccpt[j]++;
        }
      }

    }

    u = gsl_rng_uniform(gsl_r);
    if(u<fswap)
    {
      ibeta = gsl_rng_uniform_int(gsl_r, nbeta-1);
      ubeta = gsl_rng_uniform(gsl_r);

      var_con1[0] = exp(theta[ibeta+1][0]);
      var_con1[1] = exp(theta[ibeta+1][1]);
      var_con1[2] = 1.0;
      if(parset.flag_line == 1)
      {
        var_hb1[0] = exp(theta[ibeta+1][2]);
        var_hb1[1] = exp(theta[ibeta+1][3]);
        var_hb1[2] = 1.0;
      }
      set_scale(theta[ibeta+1]);
      scale_flux();
      prob1 = prob_variability_semi_beta(var_con1, var_hb1, beta[ibeta]);

      var_con2[0] = exp(theta[ibeta][0]);
      var_con2[1] = exp(theta[ibeta][1]);
      var_con2[2] = 1.0;
      if(parset.flag_line == 1)
      {
        var_hb2[0] = exp(theta[ibeta][2]);
        var_hb2[1] = exp(theta[ibeta][3]);
        var_hb2[2] = 1.0;
      }
      set_scale(theta[ibeta]);
      scale_flux();
      prob2 = prob_variability_semi_beta(var_con2, var_hb2, beta[ibeta+1]);
     
      ratio = ( prob1 + prob2 ) -( prob[ibeta] + prob[ibeta+1] );
      
      if(log(ubeta) < ratio)
      {
        memcpy(swap, theta[ibeta], sizeof(theta[0]));
        memcpy(theta[ibeta], theta[ibeta+1], sizeof(theta[0]));
        memcpy(theta[ibeta+1], swap, sizeof(theta[0]));

        prob[ibeta] = prob1;
        prob[ibeta+1] = prob2;

        printf("%d and %d swap!\n", ibeta, ibeta+1);
      }
    }
    
    for(j=0; j<nbeta; j++)
      for(i=0; i<ntheta; i++)
        theta_mcmc[j][i*n_cov_update + istep%n_cov_update] = theta[j][i];

    if(istep >= nbuilt)
    { 
      fprintf(fmcmc_out, "%d", istep);
      for(i=0;i<ntheta;i++)
      {
        if(i<ndrw*2+ncode-1)
          fprintf(fmcmc_out,"\t%f", theta[nbeta-1][i]);
        else
          fprintf(fmcmc_out,"\t%f", theta[nbeta-1][i]*flux_mean_cont); // acount for flux_mean_cont
      }
      fprintf(fmcmc_out, "\n");
    }

    if(istep%1000==0)printf("step %d\n", istep);

    istep++;

    if(istep%n_cov_update == 0)
    {
      for(j=0; j<nbeta; j++)
      {
        get_cov_matrix_diag(theta_mcmc[j], n_cov_update, ntheta);
//      display_mat(cov_matrix, ntheta, ntheta);
        //for(i=0; i<ntheta; i++)printf("%f\n", cov_matrix[i*ntheta+i]);
        memcpy(Pmatrix[j], cov_matrix, ntheta*ntheta*sizeof(double));
        Chol_decomp_U(Pmatrix[j], ntheta, &info);
      }
    }
  }
  
}

void set_mcmc_param(double *theta, double *stepsize, double (*range)[2], int ntheta)
{
  int i, j;

  i=0;
  range[i][0] = log(0.001);
  range[i++][1] = log(10.0);
  range[i][0] = log(cadence);
  range[i++][1] = log(date_span_cont);

  if(parset.flag_line == 1)
  {
    range[i][0] = log(0.001);
    range[i++][1] = log(10.0);
    range[i][0] = log(cadence);
    range[i++][1] = log(date_span_line);
  }

  for(j=0; j<ncode-1; j++)
  {
    range[i][0] = parset.scale_range_low;
    range[i++][1] = parset.scale_range_up;
  }
  for(j=0; j<ncode-1; j++)
  {
    range[i][0] = parset.shift_range_low;
    range[i++][1] = parset.shift_range_up;
  }

  if(i!=ntheta)
  {
    printf("Wrong in set_mcmc_param()!\n");
    exit(-1);
  }

  i=0;
  theta[i] = (range[i][0] + range[i][1])*0.5;   // log sigma
  i++;
  theta[i] = (range[i][0] + range[i][1])*0.5;   // log tau
  i++;
  
  if(parset.flag_line == 1)
  {
    theta[i] = (range[i][0] + range[i][1])*0.5;   // log sigma
    i++;
    theta[i] = (range[i][0] + range[i][1])*0.5;   // log tau
    i++;
  }
  
  for(j=0; j<ncode-1; j++)
  {
    theta[i++] = sqrt(parset.scale_range_low * parset.scale_range_up);
  }
  for(j=0; j<ncode-1; j++)
  {
    theta[i++] = 0.5*(parset.shift_range_low + parset.shift_range_up);
  }
  
  if(i!=ntheta)
  {
    printf("Wrong in set_mcmc_param()!\n");
    exit(-1);
  }


  i=0;
  stepsize[i++] = 0.1;
  stepsize[i++] = 0.1;

  if(parset.flag_line == 1)
  {
    stepsize[i++] = 0.1;
    stepsize[i++] = 0.1;
  }

  for(j=0; j<ncode-1; j++)
  {
    stepsize[i++] = 0.01;
  }
  for(j=0; j<ncode-1; j++)
  {
    stepsize[i++] = 0.01;
  }
  
  if(i!=ntheta)
  {
    printf("Wrong in set_mcmc_param()!\n");
    exit(-1);
  }
}

void init_cov_matrix(double *stepsize, int ntheta)
{
  int i, j;
  for(i=0; i<ntheta; i++)
  {
    cov_matrix[i*ntheta + i] = stepsize[i] * stepsize[i];
    for(j=0; j<i; j++)
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = 0.0;
  }
}

void scale_flux()
{
  int i, idx;
  for(i=0; i<nd_cont; i++)
  {
    idx = code_idx_cont[i];
    Fcon[i] = optflux[i][0] * ps_scale[idx] - es_scale[idx];
    Fcon_err[i] = optflux[i][1] * ps_scale[idx];
  }

  if(parset.flag_line == 1)
  {
    for(i=0; i<nd_line; i++)
    {
      idx = code_idx_line[i];
      Fhb[i] = hbb[i][0] * ps_scale[idx];
      Fhb_err[i] = hbb[i][1] * ps_scale[idx];
    }
  }
}

void scale_flux_err()
{
  int i, idx;
  for(i=0; i<nd_cont; i++)
  {
    idx = code_idx_cont[i];
    Fcon[i] = optflux[i][0] * ps_scale[idx] - es_scale[idx];
    //Fcon_err[i] = optflux[i][1] * ps_scale[idx];
    Fcon_err[i] = sqrt(pow(optflux[i][1] * ps_scale[idx], 2.0)
                      +pow(optflux[i][0] * ps_scale_err[idx], 2.0)
                      +pow(es_scale_err[idx], 2.0)
                      -2.0*optflux[i][0]*pe_scale_covar[idx]);
  }

  if(parset.flag_line == 1)
  {
    for(i=0; i<nd_line; i++)
    {
      idx = code_idx_line[i];
      Fhb[i] = hbb[i][0] * ps_scale[idx];
      Fhb_err[i] = sqrt(pow(hbb[i][1] * ps_scale[idx], 2) 
                      + pow(hbb[i][0] * ps_scale_err[idx], 2.0) );
    }
  }
}

void set_scale(double *theta)
{
  int i;
  for(i=1; i<ncode; i++)
  {
    ps_scale[i] = theta[ndrw*2+i-1];
    es_scale[i] = theta[ndrw*2+ncode-1+i-1];
  }

}

void set_scale_err(double *theta_var)
{
  int i;
  for(i=1; i<ncode; i++)
  {
    ps_scale_err[i] = theta_var[ndrw*2+i-1];
    es_scale_err[i] = theta_var[ndrw*2+ncode-1+i-1];
  }

}

void set_scale_covar(double **theta_covar)
{
  int i;
  for(i=1; i<ncode; i++)
  {
    pe_scale_covar[i] = theta_covar[ndrw*2+i-1][ndrw*2+ncode-1+i-1];
  }
}

void mcmc_memory_init()
{
  int i;
 
  Smat = array_malloc(nd_max*nd_max);
  Nmat = array_malloc(nd_max*nd_max);
  N0mat = array_malloc(nd_max*nd_max);
  Qmat = array_malloc(nd_max*nd_max);
  ISmat = array_malloc(nd_max*nd_max);
  INmat = array_malloc(nd_max*nd_max);
  IQmat = array_malloc(nd_max*nd_max);
  Cmat = array_malloc(nd_max*nd_max);
  ICmat = array_malloc(nd_max*nd_max);
    
  Fcon = array_malloc(nd_max);
  Fcon_err = array_malloc(nd_max);
  Fhb = array_malloc(nd_max);
  Fhb_err = array_malloc(nd_max);
  

  cov_matrix = array_malloc(ntheta_max * ntheta_max);

  var_con_best = array_malloc(3*sizeof(double));
  var_con_best_err = array_malloc(3*sizeof(double));
  var_hb_best = array_malloc(3*sizeof(double));
  var_hb_best_err = array_malloc(3*sizeof(double));

  workspace = array_malloc(10*nd_max);

//GSL
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  gsl_rng_set(gsl_r, time(NULL)); 
}

void mcmc_memory_free()
{
  free(Smat);
  free(Nmat);
  free(N0mat);
  free(Qmat);
  free(ISmat);
  free(INmat);
  free(IQmat);
  free(Cmat);
  free(ICmat);

  free(Fcon);
  free(Fcon_err);
  free(Fhb);
  free(Fhb_err);

  free(var_con_best);
  free(var_con_best_err);
  free(var_hb_best);
  free(var_hb_best_err);

  free(cov_matrix);
  free(workspace);
}

/**
 * probability for given variability model parameters.
 */
double prob_variability(double *var_con, double *var_hb)
{
  double prob, prob1=0.0, prob2=0.0, lambda, ave_con, lndet, sigma, tau, alpha;
  double lndet_n0, lndet_n, prior_phi;
  double * ybuf, * Larr;
  int i, info;

  Larr = workspace;
  ybuf = Larr + nd_max;

  tau = var_con[1];
  sigma = var_con[0] * sqrt(tau);
  alpha = var_con[2];

  set_covar_mat_con(sigma, tau, alpha);
  
  for(i=0;i<nd_cont*nd_cont; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }
  
  memcpy(ICmat, Cmat, nd_cont*nd_cont*sizeof(double));
  inverse_mat(ICmat, nd_cont, &info);
  
  /* the best estimate for average */
  for(i=0;i<nd_cont;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nd_cont, ybuf);
  lambda = cblas_ddot(nd_cont, Larr, 1, ybuf, 1);
  multiply_matvec(ICmat, Fcon, nd_cont, ybuf);
  ave_con = cblas_ddot(nd_cont, Larr, 1, ybuf, 1);
  ave_con /=lambda;

/* get the probability */
  for(i=0;i<nd_cont;i++)Larr[i] = Fcon[i] - ave_con;
  multiply_matvec(ICmat, Larr, nd_cont, ybuf);
  prob1 = -0.5 * cblas_ddot(nd_cont, Larr, 1, ybuf, 1);

  lndet = lndet_mat(Cmat, nd_cont, &info);
//  lndet_n = lndet_mat(Nmat, nd_cont, &info);
//  lndet_n0 = lndet_mat(N0mat, nd_cont, &info);

  lndet_n = lndet_n0 = 0.0;
  for(i=0; i<nd_cont; i++)
  {
    lndet_n += 2.0*log(Fcon_err[i]);
    lndet_n0 += 2.0*log(optflux[i][1]);
  }

  prob1 = prob1 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);

  if(parset.flag_line == 1)
  {
    tau = var_hb[1];
    sigma = var_hb[0]*sqrt(tau);
    alpha = var_hb[2];

    set_covar_mat_hb(sigma, tau, alpha);
  
    for(i=0;i<nd_line*nd_line; i++)
    {
      Cmat[i] = Smat[i] + Nmat[i];
    }
  
    memcpy(ICmat, Cmat, nd_line*nd_line*sizeof(double));
    inverse_mat(ICmat, nd_line, &info);
  
  /* the best estimate for average */
    for(i=0;i<nd_line;i++)Larr[i]=1.0;
    multiply_matvec(ICmat, Larr, nd_line, ybuf);
    lambda = cblas_ddot(nd_line, Larr, 1, ybuf, 1);
    multiply_matvec(ICmat, Fhb, nd_line, ybuf);
    ave_con = cblas_ddot(nd_line, Larr, 1, ybuf, 1);
    ave_con /=lambda;

/* get the probability */
    for(i=0;i<nd_line;i++)Larr[i] = Fhb[i] - ave_con;
    multiply_matvec(ICmat, Larr, nd_line, ybuf);
    prob2 = -0.5 * cblas_ddot(nd_line, Larr, 1, ybuf, 1);

    lndet = lndet_mat(Cmat, nd_line, &info);
//    lndet_n = lndet_mat(Nmat, nd_line, &info);
//    lndet_n0 = lndet_mat(N0mat, nd_line, &info);
    lndet_n = lndet_n0 = 0.0;
    for(i=0; i<nd_line; i++)
    {
      lndet_n += 2.0*log(Fhb_err[i]);
      lndet_n0 += 2.0*log(hbb[i][1]);
    }
    prob2 = prob2 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
//  prob = prob - 0.5*log(lndet);
  }

  prior_phi = 1.0;
  for(i=1; i<ncode; i++)
  {
    prior_phi *= ps_scale[i]; 
  }
  
  prob = (prob1 + prob2 - log(prior_phi));
  
  return prob;
}

/**
 * probability for given variability model parameters.
 */
double prob_variability_beta(double *var_con, double *var_hb, double beta)
{
  double prob, prob1=0.0, prob2=0.0, lambda, ave_con, lndet, sigma, tau, alpha;
  double lndet_n, lndet_n0, prior_phi;
  double * ybuf, * Larr;
  int i, info;

  Larr = workspace;
  ybuf = Larr + nd_cont;

  tau = var_con[1];
  sigma = var_con[0] * sqrt(tau);
  alpha = var_con[2];

  set_covar_mat_con(sigma, tau, alpha);
  
  for(i=0;i<nd_cont*nd_cont; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }
  
  memcpy(ICmat, Cmat, nd_cont*nd_cont*sizeof(double));
  inverse_mat(ICmat, nd_cont, &info);
  
  /* the best estimate for average */
  for(i=0;i<nd_cont;i++)Larr[i]=1.0;
  multiply_matvec(ICmat, Larr, nd_cont, ybuf);
  lambda = cblas_ddot(nd_cont, Larr, 1, ybuf, 1);
  multiply_matvec(ICmat, Fcon, nd_cont, ybuf);
  ave_con = cblas_ddot(nd_cont, Larr, 1, ybuf, 1);
  ave_con /=lambda;
//  printf("%f\n", ave_con);

/* get the probability */
  for(i=0;i<nd_cont;i++)Larr[i] = Fcon[i] - ave_con;
  multiply_matvec(ICmat, Larr, nd_cont, ybuf);
  prob1 = -0.5 * cblas_ddot(nd_cont, Larr, 1, ybuf, 1);

  lndet = lndet_mat(Cmat, nd_cont, &info);
  //lndet_n = lndet_mat(Nmat, nd_cont, &info);
  //lndet_n0 = lndet_mat(N0mat, nd_cont, &info);
  lndet_n = lndet_n0 = 0.0;
  for(i=0; i<nd_cont; i++)
  {
    lndet_n += 2.0*log(Fcon_err[i]);
    lndet_n0 += 2.0*log(optflux[i][1]);
  }
  prob1 = prob1 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  
  if(parset.flag_line == 1)
  {
    tau = var_hb[1];
    sigma = var_hb[0]*sqrt(tau);
    alpha = var_hb[2];

    set_covar_mat_hb(sigma, tau, alpha);
  
    for(i=0;i<nd_line*nd_line; i++)
    {
      Cmat[i] = Smat[i] + Nmat[i];
    }
  
    memcpy(ICmat, Cmat, nd_line*nd_line*sizeof(double));
    inverse_mat(ICmat, nd_line, &info);
  
  /* the best estimate for average */
    for(i=0;i<nd_line;i++)Larr[i]=1.0;
    multiply_matvec(ICmat, Larr, nd_line, ybuf);
    lambda = cblas_ddot(nd_line, Larr, 1, ybuf, 1);
    multiply_matvec(ICmat, Fhb, nd_line, ybuf);
    ave_con = cblas_ddot(nd_line, Larr, 1, ybuf, 1);
    ave_con /=lambda;
//  printf("%f\n", ave_con);

/* get the probability */
    for(i=0;i<nd_line;i++)Larr[i] = Fhb[i] - ave_con;
    multiply_matvec(ICmat, Larr, nd_line, ybuf);
    prob2 = -0.5 * cblas_ddot(nd_line, Larr, 1, ybuf, 1);

    lndet = lndet_mat(Cmat, nd_line, &info);
//    lndet_n = lndet_mat(Nmat, nd_line, &info);
//    lndet_n0 = lndet_mat(N0mat, nd_line, &info);
    lndet_n = lndet_n0 = 0.0;
    for(i=0; i<nd_line; i++)
    {
      lndet_n += 2.0*log(Fhb_err[i]);
      lndet_n0 += 2.0*log(hbb[i][1]);
    }

    prob2 = prob2 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  }
  
  prior_phi = 1.0;
  for(i=1; i<ncode; i++)
  {
    prior_phi *= ps_scale[i]; 
  }
  
  prob = beta*(prob1 + prob2 - log(prior_phi));
  
  return prob;
}

/**
 * probability for given variability model parameters.
 */
double prob_variability_semi(double *var_con, double *var_hb)
{
  double prob, prob1=0.0, prob2=0.0, lambda, ave_con, lndet, sigma, sigma2, tau, alpha;
  double lndet_n, lndet_n0, prior_phi;
  double * ybuf, * Larr, *W, *D, *phi, *Cq, *Lbuf, *yq;
  int i, nq, info;

  nq = 1;
  Larr = workspace;
  Lbuf = Larr + nd_cont*nq;
  ybuf = Lbuf + nd_cont*nq;
  W = ybuf + nd_cont;
  D = W + nd_cont;
  phi = D + nd_cont;
  Cq = phi + nd_cont;
  yq = Cq + nq*nq;

  tau = var_con[1];
  sigma = var_con[0] * sqrt(tau);
  sigma2 = sigma*sigma;
  alpha = var_con[2];
  
  for(i=0;i<nd_cont;i++)
    Larr[i]=1.0;

  compute_semiseparable_drw(date_cont, nd_cont, sigma2, 1.0/tau, Fcon_err, 0.0, W, D, phi);
  lndet = 0.0;
  for(i=0; i<nd_cont; i++)
    lndet += log(D[i]);


  /* calculate L^T*C^-1*L */
  multiply_mat_semiseparable_drw(Larr, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, nd_cont);

  /* calculate L^T*C^-1*y */
  multiply_matvec_semiseparable_drw(Fcon, W, D, phi, nd_cont, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, nd_cont);
  
  lambda = Cq[0];
  ave_con = yq[0]/Cq[0];

/* get the probability */
  for(i=0;i<nd_cont;i++)
  {
    ybuf[i] = Fcon[i] - ave_con;
  }
  multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_cont, sigma2, Lbuf);
  prob1 = -0.5 * cblas_ddot(nd_cont, ybuf, 1, Lbuf, 1);

  lndet_n = lndet_n0 = 0.0;
  for(i=0; i<nd_cont; i++)
  {
    lndet_n += 2.0*log(Fcon_err[i]);
    lndet_n0 += 2.0*log(optflux[i][1]);
  }
  prob1 = prob1 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  
  if(parset.flag_line == 1)
  {
    Larr = workspace;
    Lbuf = Larr + nd_line*nq;
    ybuf = Lbuf + nd_line*nq;
    W = ybuf + nd_line;
    D = W + nd_line;
    phi = D + nd_line;
    Cq = phi + nd_line;
    yq = Cq + nq*nq;

    tau = var_hb[1];
    sigma = var_hb[0]*sqrt(tau);
    alpha = var_hb[2];

    for(i=0;i<nd_line;i++)
      Larr[i]=1.0;

    compute_semiseparable_drw(date_line, nd_line, sigma2, 1.0/tau, Fhb_err, 0.0, W, D, phi);
    lndet = 0.0;
    for(i=0; i<nd_line; i++)
      lndet += log(D[i]);

    /* calculate L^T*C^-1*L */
    multiply_mat_semiseparable_drw(Larr, W, D, phi, nd_line, nq, sigma2, Lbuf);
    multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, nd_line);

    /* calculate L^T*C^-1*y */
    multiply_matvec_semiseparable_drw(Fhb, W, D, phi, nd_line, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, nd_line);
  
    lambda = Cq[0];
    ave_con = yq[0]/Cq[0];

    for(i=0;i<nd_line;i++)
    {
      ybuf[i] = Fhb[i] - ave_con;
    }
    multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_line, sigma2, Lbuf);
    prob2 = -0.5 * cblas_ddot(nd_line, ybuf, 1, Lbuf, 1);

    lndet_n = lndet_n0 = 0.0;
    for(i=0; i<nd_line; i++)
    {
      lndet_n += 2.0*log(Fhb_err[i]);
      lndet_n0 += 2.0*log(hbb[i][1]);
    }

    prob2 = prob2 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  }
  
  prior_phi = 1.0;
  for(i=1; i<ncode; i++)
  {
    prior_phi *= ps_scale[i]; 
  }
  
  prob = prob1 + prob2 - log(prior_phi);
  
  return prob;
}

/**
 * probability for given variability model parameters.
 */
double prob_variability_semi_beta(double *var_con, double *var_hb, double beta)
{
  double prob, prob1=0.0, prob2=0.0, lambda, ave_con, lndet, sigma, sigma2, tau, alpha;
  double lndet_n, lndet_n0, prior_phi;
  double * ybuf, * Larr, *W, *D, *phi, *Cq, *Lbuf, *yq;
  int i, nq, info;

  nq = 1;
  Larr = workspace;
  Lbuf = Larr + nd_cont*nq;
  ybuf = Lbuf + nd_cont*nq;
  W = ybuf + nd_cont;
  D = W + nd_cont;
  phi = D + nd_cont;
  Cq = phi + nd_cont;
  yq = Cq + nq*nq;

  tau = var_con[1];
  sigma = var_con[0] * sqrt(tau);
  sigma2 = sigma*sigma;
  alpha = var_con[2];
  
  for(i=0;i<nd_cont;i++)
    Larr[i]=1.0;

  compute_semiseparable_drw(date_cont, nd_cont, sigma2, 1.0/tau, Fcon_err, 0.0, W, D, phi);
  lndet = 0.0;
  for(i=0; i<nd_cont; i++)
    lndet += log(D[i]);


  /* calculate L^T*C^-1*L */
  multiply_mat_semiseparable_drw(Larr, W, D, phi, nd_cont, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, nd_cont);

  /* calculate L^T*C^-1*y */
  multiply_matvec_semiseparable_drw(Fcon, W, D, phi, nd_cont, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, nd_cont);
  
  lambda = Cq[0];
  ave_con = yq[0]/Cq[0];

/* get the probability */
  for(i=0;i<nd_cont;i++)
  {
    ybuf[i] = Fcon[i] - ave_con;
  }
  multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_cont, sigma2, Lbuf);
  prob1 = -0.5 * cblas_ddot(nd_cont, ybuf, 1, Lbuf, 1);

  lndet_n = lndet_n0 = 0.0;
  for(i=0; i<nd_cont; i++)
  {
    lndet_n += 2.0*log(Fcon_err[i]);
    lndet_n0 += 2.0*log(optflux[i][1]);
  }
  prob1 = prob1 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  
  if(parset.flag_line == 1)
  {
    Larr = workspace;
    Lbuf = Larr + nd_line*nq;
    ybuf = Lbuf + nd_line*nq;
    W = ybuf + nd_line;
    D = W + nd_line;
    phi = D + nd_line;
    Cq = phi + nd_line;
    yq = Cq + nq*nq;
    
    tau = var_hb[1];
    sigma = var_hb[0]*sqrt(tau);
    alpha = var_hb[2];

    for(i=0;i<nd_line;i++)
      Larr[i]=1.0;

    compute_semiseparable_drw(date_line, nd_line, sigma2, 1.0/tau, Fhb_err, 0.0, W, D, phi);
    lndet = 0.0;
    for(i=0; i<nd_line; i++)
      lndet += log(D[i]);

    /* calculate L^T*C^-1*L */
    multiply_mat_semiseparable_drw(Larr, W, D, phi, nd_line, nq, sigma2, Lbuf);
    multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, nd_line);

    /* calculate L^T*C^-1*y */
    multiply_matvec_semiseparable_drw(Fhb, W, D, phi, nd_line, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, nd_line);
  
    lambda = Cq[0];
    ave_con = yq[0]/Cq[0];

    for(i=0;i<nd_line;i++)
    {
      ybuf[i] = Fhb[i] - ave_con;
    }
    multiply_matvec_semiseparable_drw(ybuf, W, D, phi, nd_line, sigma2, Lbuf);
    prob2 = -0.5 * cblas_ddot(nd_line, ybuf, 1, Lbuf, 1);

    lndet_n = lndet_n0 = 0.0;
    for(i=0; i<nd_line; i++)
    {
      lndet_n += 2.0*log(Fhb_err[i]);
      lndet_n0 += 2.0*log(hbb[i][1]);
    }

    prob2 = prob2 - 0.5*lndet - 0.5*log(lambda) + 0.5 * (lndet_n - lndet_n0);
  }
  
  prior_phi = 1.0;
  for(i=1; i<ncode; i++)
  {
    prior_phi *= ps_scale[i]; 
  }
  
  prob = beta*(prob1 + prob2 - log(prior_phi));
  
  return prob;
}

void set_covar_mat_con(double sigma, double tau, double alpha)
{
  double t1, t2, nerr;
  int i, j;

  for(i=0; i<nd_cont; i++)
  {
    t1 = date_cont[i];
    for(j=0; j<=i; j++)
    {
      t2 = date_cont[j];
      Smat[i*nd_cont+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha)); 
      Smat[j*nd_cont+i] = Smat[i*nd_cont+j];
      
      Nmat[i*nd_cont+j] = Nmat[j*nd_cont+i] = 0.0; 
      
      N0mat[i*nd_cont+j] = N0mat[j*nd_cont+i] = 0.0; 
    }
    nerr = Fcon_err[i];
    Nmat[i*nd_cont+i] = nerr*nerr;
    
    N0mat[i*nd_cont + i] = optflux[i][1] * optflux[i][1];
  }

  return;
}

void set_covar_mat_hb(double sigma, double tau, double alpha)
{
  double t1, t2, nerr;
  int i, j;

  for(i=0; i<nd_line; i++)
  {
    t1 = date_line[i];
    for(j=0; j<=i; j++)
    {
      t2 = date_line[j];
      Smat[i*nd_line+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha)); 
      Smat[j*nd_line+i] = Smat[i*nd_line+j];
      
      Nmat[i*nd_line+j] = Nmat[j*nd_line+i] = 0.0; 
      N0mat[i*nd_line+j] = N0mat[j*nd_line+i] = 0.0; 
    }
    nerr = Fhb_err[i];
    Nmat[i*nd_line+i] = nerr*nerr;
    N0mat[i*nd_line + i] = hbb[i][1] * hbb[i][1];
  }

  return;
}

void mcmc_stats()
{
  FILE *fmcmc, *fout;
  int i, j, ip, nstep, ntheta;
  double **theta, *theta_mean, *theta_var, **theta_covar;
  const int nh=100;    
  gsl_histogram **hd;
  
  ntheta = ndrw*2 + 2*(ncode-1);
  hd = malloc(ntheta*sizeof(gsl_histogram *));
  for(i=0; i<ntheta; i++)
  {
    hd[i] = gsl_histogram_alloc(nh);
  }
  
  theta = matrix_malloc(ntheta, n_mcmc-nbuilt);
  theta_mean = array_malloc(ntheta);
  theta_var = array_malloc(ntheta);
  theta_covar = matrix_malloc(ntheta, ntheta);
  
  fmcmc = fopen("mcmc.txt", "r");
  nstep = 0;
  while(!feof(fmcmc) && nstep<n_mcmc-nbuilt-1)
  {
    fscanf(fmcmc,"%d", &ip);
    for(i=0; i<ntheta; i++)
    {
      fscanf(fmcmc,"%lf", &theta[i][nstep]);
      if(i >= ndrw*2+ncode-1)
        theta[i][nstep] /= flux_mean_cont; // acount for flux_mean_cont;
    }
    fscanf(fmcmc,"\n");

    nstep++;
    if(nstep>n_mcmc-nbuilt)
    {
      fprintf(stderr, "nstep > n_mcmc");
      exit(-1);
    }
  }
  
  for(i=0; i<ntheta; i++)
  {
    theta_mean[i] = gsl_stats_mean(&theta[i][0], 1, nstep);
    theta_var[i] = gsl_stats_sd(&theta[i][0], 1, nstep);
    for(j=0; j<i; j++)
    {
      theta_covar[i][j] = theta_covar[j][i] = gsl_stats_covariance(&theta[i][0], 1, &theta[j][0], 1, nstep);
    }
    theta_covar[i][i] = theta_var[i] * theta_var[i];
  }
  
  for(i=0; i<ntheta; i++)
  {
    printf("%f %f\n", theta_mean[i], theta_var[i]);
  }
  set_scale(theta_mean);
  set_scale_err(theta_var);
  set_scale_covar(theta_covar);
  scale_flux_err();
  
  fout = fopen("factor.txt", "w");
  printf("factor:\n");
  for(i=0; i<ncode; i++)
  {
    printf("%s\t%f\t%f\t%e\n", code[i], ps_scale[i], es_scale[i]*flux_mean_cont, pe_scale_covar[i]*flux_mean_cont);
    fprintf(fout, "%f\t%f\t%f\t%f\t%e\t%s\n", ps_scale[i], ps_scale_err[i], 
      es_scale[i]*flux_mean_cont, es_scale_err[i]*flux_mean_cont,  pe_scale_covar[i]*flux_mean_cont, code[i]);
  }
  
  for(i=0; i<2; i++)
  {
    var_con_best[i] = exp(theta_mean[i]);
    var_con_best_err[i] = theta_var[i] * var_con_best[i];
  }
  var_con_best[2] = 1.0;
  var_hb_best[2] = 1.0;

  if(parset.flag_line == 1)
  {
    for(i=0; i<2; i++)
    {
      var_hb_best[i] = exp(theta_mean[i+2]);
      var_hb_best_err[i] = theta_var[i+2] * var_hb_best[i];
    }
    var_con_best_err[2] = 0.0;
    var_hb_best_err[2] = 0.0;  
  }

  //printf("%f\t%f\t%f\n", var_con_best[0], var_con_best[1], var_con_best[2]);
  //printf("%f\t%f\t%f\n", var_hb_best[0], var_hb_best[1], var_hb_best[2]);

  for(i=0; i<ntheta; i++)
    free(theta[i]);
  free(theta);

  free(theta_mean);
  free(theta_var);

  for(i=0; i<ntheta; i++)
    free(theta_covar[i]);
  free(theta_covar);

  fclose(fmcmc);
  fclose(fout);
  return;
}

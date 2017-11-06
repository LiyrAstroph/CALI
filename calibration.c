#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int calibrate()
{

  n_mcmc = 100000;
  nbuilt = 50000;
  n_recon = 500;
  mcmc_memory_init();
  reconstruct_int();
  //mcmc_pt();
  mcmc_stats();

  printf("%f\t%f\t%f\n", var_con_best[0], var_con_best[1], var_con_best[2]);
  printf("%f\t%f\t%f\n", var_hb_best[0], var_hb_best[1], var_hb_best[2]);

  reconstruct_con(var_con_best);
  reconstruct_hb(var_hb_best);
}


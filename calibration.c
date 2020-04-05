#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int calibrate()
{

  n_mcmc = parset.n_mcmc;
  nbuilt = parset.n_builtin;
  n_recon = 500;
  mcmc_memory_init();
  reconstruct_int();
  mcmc_pt();
  mcmc_stats();

  output_optical("cont.txt");
  output_hb("line.txt");

  printf("Cont DRW:ln(sigma)=%f\tln(tau)=%f\n", var_con_best[0], var_con_best[1]);
  reconstruct_con(var_con_best);

  if(parset.flag_line == 1)
  {
    printf("Line DRW:ln(sigma)=%f\tln(tau)=%f\n", var_hb_best[0], var_hb_best[1]);
    reconstruct_hb(var_hb_best);
  }
  reconstruct_end();

  mcmc_memory_free();

  return EXIT_SUCCESS;
}


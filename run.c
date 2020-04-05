/*
 * code for relative scaling and inter-calibration of spectra
 * used in the paper Li, Y.-R. et al., 2013 
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "proto.h"
#include "allvars.h"

int run(char *file_param)
{
  strcpy(parset.file_param, file_param);

  read_parset();

  memory_malloc();
  read_dataset();
  
  calibrate();
   
  memory_free();

  return EXIT_SUCCESS;
}
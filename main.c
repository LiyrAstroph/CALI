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
#include <ctype.h>

#include "proto.h"
#include "allvars.h"


//start of the program
int main(int argc, char *argv[])
{
  if(argc <2)
  {
    printf("# No parameter file is specified.\n");
    exit(0);
  }
  strcpy(parset.file_param, argv[1]);

  read_parset();

  memory_malloc();
  read_dataset();
  
  calibrate();
   
  memory_free();
  return 0;
}
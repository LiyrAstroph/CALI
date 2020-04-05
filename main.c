/*
 * code for relative scaling and inter-calibration of spectra
 * used in the paper Li, Y.-R. et al., 2013 
 *
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 */

#include <stdio.h>
#include <stdlib.h>

int run(char *file_param);

//start of the program
int main(int argc, char *argv[])
{
  if(argc <2)
  {
    printf("# No parameter file is specified.\n");
    exit(0);
  }
  run(argv[1]);
  return EXIT_SUCCESS;
}
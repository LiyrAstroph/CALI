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
  memory_malloc();
  read_dataset("ngc5548_year1.txt");

/*  get_file_list("contifit/year1_sort.txt");
  read_dateset_code("group_year1.txt");
  get_date_from_filename();
  
  get_cont_fitting_par("result/year1_sort.txt.conti");
  get_line_fitting_par("result/year1_sort.txt.hbo3_gh4");
  get_o3flux();
  get_o3center();
  scale_flux_by_OIII();*/
  
  calibrate();
  
  output_optical("opt.txt");
  output_hb("hb.txt");
 
//  scaling();
  return 0;
}

int scaling()
{
  int i, j, nt, eof;
  char fname[100], *pstr, str[100];
  const int nd_max=50000;
  double *wave, *flux, scale, shift;
  FILE *fp, *fs, *flist;
  
  wave = malloc(nd_max * sizeof(double));
  flux = malloc(nd_max * sizeof(double));
  
  flist = fopen("list_scaled.txt", "w");

  for(i=0; i<nlist; i++)
  {
    strcpy(fname, "orig/");
    strcat(fname, list[i]);
    printf("%s\n", fname);

    scale = o3flux_std / o3flux[i];
    shift = 5006.84 - o3center[i];

    fp = fopen(fname, "r");
    for(j=0; j<nd_max; j++)
    {
      eof = fscanf(fp, "%lf %lf", &wave[j], &flux[j]);
      if(eof<=0)
        break;
    }
    nt = j;
    fclose(fp);

    for(j=0; j<nt; j++)
    {
      wave[j] += shift;
      flux[j] *= scale;
    }

    strcpy(fname, "orig/");
    strcat(fname, list[i]);
    pstr=strrchr(fname, '.');
    *pstr='\0';
    strcat(fname, "_scaled.txt");

    printf("%s\n", fname);
    fs = fopen(fname, "w");

    for(j=0; j<nt; j++)
    {
      fprintf(fs, "%f\t%f\n", wave[j], flux[j]);
    } 

    fclose(fs);
    
    strcpy(str, list[i]);
    pstr=strrchr(str, '.');
    *pstr='\0';
    strcat(str, "_scaled.txt");
    fprintf(flist, "%s     %8.6f     %d     %s\n", str, 0.017175,    0,    "1d0");
  }
  
  fclose(flist);
}


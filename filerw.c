#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

int memory_malloc()
{
  int i;

  list = malloc(nd_max*sizeof(char *));
  for(i=0; i<nd_max; i++)
  {
    list[i] = malloc(100*sizeof(char));
  }

// emission line  
  hbb = malloc(nd_max * sizeof(double *));
  o3b = malloc(nd_max * sizeof(double *));
  o3n = malloc(nd_max * sizeof(double *));
  hbn = malloc(nd_max * sizeof(double *));
  chi2 = malloc(nd_max * sizeof(double));

  for(i=0; i<nd_max; i++)
  {
    hbb[i] = malloc(10*sizeof(double));
    o3b[i] = malloc(6*sizeof(double));
    o3n[i] = malloc(6*sizeof(double));
    hbn[i] = malloc(4*sizeof(double));
  }

  o3flux = malloc(nd_max * sizeof(double));
  o3center = malloc(nd_max * sizeof(double));

// continuum
  optflux = malloc(nd_max * sizeof(double *));
  optslope = malloc(nd_max * sizeof(double *));
  chi2opt = malloc(nd_max * sizeof(double));
  
  for(i=0; i<nd_max; i++)
  {
    optflux[i] = malloc(2*sizeof(double));
    optslope[i] = malloc(2*sizeof(double));
  }


  date_cont = malloc(nd_max * sizeof(double));
  date_line = malloc(nd_max * sizeof(double));
  code = malloc(ncode_max * sizeof(char *));
  for(i=0; i<ncode_max; i++)
  {
    code[i] = malloc(100 * sizeof(char));
  }
  obs_num_cont = malloc(ncode_max * sizeof(int));
  obs_num_line = malloc(ncode_max * sizeof(int));
  code_idx_cont = malloc(nd_max * sizeof(int));
  code_idx_line = malloc(nd_max * sizeof(int));
  
  ps_scale = malloc(ncode_max * sizeof(double));
  es_scale = malloc(ncode_max * sizeof(double));

  ps_scale_err = malloc(ncode_max * sizeof(double));
  es_scale_err = malloc(ncode_max * sizeof(double));
  
  for(i=0; i<ncode_max; i++)
  {
    ps_scale[i] = 1.0;
    es_scale[i] = 0.0;

    ps_scale_err[i] = 0.0;
    es_scale_err[i] = 0.0;
  }
  
}

int output_optical(char *fname)
{
  FILE *fout;
  int i, idx;

  fout=fopen(fname,"w");
  if(fout==NULL)
  {
    printf("cannot open file %s\n", fname);
    exit(-1);
  }

  for(i=0; i<nd_cont; i++)
  {
    idx = code_idx_cont[i];
    fprintf(fout, "%f\t%f\t%f\t%s\n", date_cont[i], optflux[i][0] * ps_scale[idx] - es_scale[idx],
     optflux[i][1]* ps_scale[idx], code[code_idx_cont[i]]);
  }

  fclose(fout);
}

int output_hb(char *fname)
{
  FILE *fout;
  int i, idx;
  
  fout=fopen(fname,"w");
  if(fout==NULL)
  {
    printf("cannot open file %s\n", fname);
    exit(-1);
  }

  for(i=0; i<nd_line; i++)
  {
    idx = code_idx_line[i];
    fprintf(fout, "%f\t%f\t%f\t%s\n", date_line[i], hbb[i][0] * ps_scale[idx], 
           hbb[i][1]* ps_scale[idx], code[code_idx_line[i]]);
  }

  fclose(fout);
}

int read_dataset()
{
  FILE *fcont, *fline;
  char str[300], strcode[100];
  int i, j, idx, ic;
  
  // read continuum 
  fcont =fopen(parset.file_cont, "r");
  if(fcont==NULL)
  {
    printf("cannot open file %s\n", parset.file_cont);
    exit(-1);
  }
  
  idx = 0;
  ic = 0;
  for(i=0; i<nd_max; i++)
  {
    fgets(str, 300, fcont);
    if(feof(fcont))
      break;
    if(sscanf(str, "# %s %d\n", code[idx], &obs_num_cont[idx]) < 2)
    {
      printf("# Wrong in reading %s.\n", parset.file_cont);
      exit(0);
    }
    
    printf("# %s %d\n", code[idx], obs_num_cont[idx]);
    for(j=0; j<obs_num_cont[idx]; j++)
    {
      fgets(str, 300, fcont);
//      printf("%s\n", str);
      if(sscanf(str, "%lf %lf %lf\n", &date_cont[ic], &optflux[ic][0], &optflux[ic][1]) < 3)
      {
        printf("# Wrong in reading %s.\n", parset.file_cont);
        exit(0);
      }
//      printf("%lf %lf %lf\n", date[ic], optflux[ic][0], hbb[ic][0]);
      code_idx_cont[ic] = idx;
      
      optflux[ic][0] *= (o3flux_std);
      optflux[ic][1] *= (o3flux_std);
      ic++;
    }
    idx++;
  }
  nd_cont = ic;
  ncode = idx;
  printf("Cont points: %d, codes: %d\n", ic, idx); 
  fclose(fcont);

  
  //read line
  if(parset.flag_line == 1)
  {
  fline =fopen(parset.file_line, "r");
  if(fline==NULL)
  {
    printf("cannot open file %s\n", parset.file_line);
    exit(-1);
  }
  
  idx = 0;
  ic = 0;
  for(i=0; i<nd_max; i++)
  {
    fgets(str, 300, fline);
    if(feof(fline))
      break;
    if(sscanf(str, "# %s %d\n", strcode, &obs_num_line[idx]) < 2)
    {
      printf("# Wrong in reading %s.\n", parset.file_line);
      exit(0);
    }
    if(strcmp(strcode, code[idx]) !=0)
    {
      printf("# code of line %s is different from cont %s\n", strcode, code[idx]);
      exit(0);
    }
    printf("# %s %d\n", code[idx], obs_num_line[idx]);
    for(j=0; j<obs_num_line[idx]; j++)
    {
      fgets(str, 300, fline);
//      printf("%s\n", str);
      if(sscanf(str, "%lf %lf %lf\n", &date_line[ic], &hbb[ic][0], &hbb[ic][1]) < 3)
      {
        printf("# Wrong in reading %s.\n", parset.file_line);
        exit(0);
      }
//      printf("%lf %lf %lf\n", date[ic], optflux[ic][0], hbb[ic][0]);
      code_idx_line[ic] = idx;

      hbb[ic][0] *= o3flux_std;
      hbb[ic][1] *= o3flux_std;

      ic++;
    }
    idx++;
  }

  if(idx != ncode)
  {
    printf("# Numbers of codes for cont and line are not same: %d %d!\n", ncode, idx);
    exit(0);
  }

  nd_line = ic;
  ncode = idx;
  printf("Line points: %d, codes: %d\n", ic, idx); 
  fclose(fline);
  }

}


/*!
 * read parameter set from parameter file.
 */
void read_parset()
{
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    nt = 0;
    strcpy(tag[nt], "FileCont");
    addr[nt] = &parset.file_cont;
    id[nt++] = STRING;

    strcpy(tag[nt], "FileLine");
    addr[nt] = &parset.file_line;
    id[nt++] = STRING;

    strcpy(tag[nt], "NMcmc");
    addr[nt] = &parset.n_mcmc;
    id[nt++] = INT;

    strcpy(tag[nt], "NBuiltin");
    addr[nt] = &parset.n_builtin;
    id[nt++] = INT;
    
    char fname[200];
    sprintf(fname, "%s", parset.file_param);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }


    /* set default configurations */
    strcpy(parset.file_cont, "");
    strcpy(parset.file_line, "");
    parset.flag_line = 0;
    parset.n_mcmc = 100000;
    parset.n_builtin = 50000;

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='#')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
        {
          j = i;
          tag[i][0] = 0;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(id[j])
        {
          case DOUBLE:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy(addr[j], buf2);
            break;
          case INT:
            *((int *)addr[j]) = (int) atof(buf2);
            break;
        }
      }
      else
      {
        fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                      parset.file_param, buf1);
        exit(0);
      }
    }
    fclose(fparam);
  
  if(parset.n_mcmc < 10000 | parset.n_builtin < 10000)
  {
    printf("# n_mcmc and n_builtin should be larger than 10000. \n");
    exit(0);
  }

  if(parset.n_mcmc <= parset.n_builtin )
  {
    printf("# n_mcmc should be larger than n_builtin.\n");
    exit(0);
  }

  if((strcmp(parset.file_line, "")==0) & (strcmp(parset.file_cont, "")==0))
  {
    fprintf(stderr, "Please specify data filename in %s with the following format:\n FileCont   xxxx\n FileLine   xxxx\n", parset.file_param);
    exit(-1);
  }
  
  if(strcmp(parset.file_line, "")!=0)
  {
    parset.flag_line = 1;
  }
  return;
}
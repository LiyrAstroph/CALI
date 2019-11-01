#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sort.h>

#include "allvars.h"
#include "proto.h"

int memory_malloc()
{
  int i;

// emission line  
  hbb = malloc(nd_max * sizeof(double *));
  for(i=0; i<nd_max; i++)
  {
    hbb[i] = malloc(2*sizeof(double));
  }

  hbb_org = malloc(nd_max * sizeof(double *));
  for(i=0; i<nd_max; i++)
  {
    hbb_org[i] = malloc(2*sizeof(double));
  }

// continuum
  optflux = malloc(nd_max * sizeof(double *));
  for(i=0; i<nd_max; i++)
  {
    optflux[i] = malloc(2*sizeof(double));
  }

  optflux_org = malloc(nd_max * sizeof(double *));
  for(i=0; i<nd_max; i++)
  {
    optflux_org[i] = malloc(2*sizeof(double));
  }


  date_cont = malloc(nd_max * sizeof(double));
  date_line = malloc(nd_max * sizeof(double));

  date_cont_org = malloc(nd_max * sizeof(double));
  date_line_org = malloc(nd_max * sizeof(double));

  perm_cont = malloc(nd_max * sizeof(int));
  perm_line = malloc(nd_max * sizeof(int));

  code = malloc(ncode_max * sizeof(char *));
  for(i=0; i<ncode_max; i++)
  {
    code[i] = malloc(100 * sizeof(char));
  }
  obs_num_cont = malloc(ncode_max * sizeof(int));
  obs_num_line = malloc(ncode_max * sizeof(int));
  code_idx_cont = malloc(nd_max * sizeof(int));
  code_idx_line = malloc(nd_max * sizeof(int));

  code_idx_cont_org = malloc(nd_max * sizeof(int));
  code_idx_line_org = malloc(nd_max * sizeof(int));

  
  ps_scale = malloc(ncode_max * sizeof(double));
  es_scale = malloc(ncode_max * sizeof(double));

  ps_scale_err = malloc(ncode_max * sizeof(double));
  es_scale_err = malloc(ncode_max * sizeof(double));
  pe_scale_covar = malloc(ncode_max * sizeof(double));

  for(i=0; i<ncode_max; i++)
  {
    ps_scale[i] = 1.0;
    es_scale[i] = 0.0;

    ps_scale_err[i] = 0.0;
    es_scale_err[i] = 0.0;

    pe_scale_covar[i] = 0.0;
  }
  
}

int memory_free()
{
  free(date_cont);
  free(date_line);
  free(date_cont_org);
  free(date_line_org);
  
  free(obs_num_cont);
  free(obs_num_line);

  free(perm_cont);
  free(perm_line);

  free(code_idx_cont);
  free(code_idx_line);
  free(code_idx_cont_org);
  free(code_idx_line_org);

  free(ps_scale);
  free(ps_scale_err);
  free(es_scale);
  free(es_scale_err);
  free(pe_scale_covar);
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
    fprintf(fout, "%f\t%f\t%f\t%s\n", date_cont[i], Fcon[i]*flux_mean_cont,
                                      Fcon_err[i]*flux_mean_cont, code[code_idx_cont[i]]);
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
    fprintf(fout, "%f\t%f\t%f\t%s\n", date_line[i], Fhb[i]*flux_mean_line, 
                                      Fhb_err[i]*flux_mean_line, code[code_idx_line[i]]);
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
      if(sscanf(str, "%lf %lf %lf\n", &date_cont_org[ic], &optflux_org[ic][0], &optflux_org[ic][1]) < 3)
      {
        printf("# Wrong in reading %s.\n", parset.file_cont);
        exit(0);
      }
//      printf("%lf %lf %lf\n", date[ic], optflux[ic][0], hbb[ic][0]);
      code_idx_cont_org[ic] = idx;
      
      ic++;
    }
    idx++;
  }
  nd_cont = ic;
  ncode = idx;
  printf("Cont points: %d, codes: %d\n", ic, idx); 
  fclose(fcont);

  flux_mean_cont = 0.0;
  idx = 0;
  for(i=0; i<obs_num_cont[idx]; i++)
    flux_mean_cont += optflux_org[i][0];
  flux_mean_cont /= obs_num_cont[idx];
  printf("Mean cont flux of code %s: %f\n", code[idx], flux_mean_cont);

  // sort over the data 
  gsl_sort_index(perm_cont, date_cont_org, 1, nd_cont);
  for(i=0; i<nd_cont; i++)
  {
    j = perm_cont[i];
    date_cont[i] = date_cont_org[j];
    code_idx_cont[i] = code_idx_cont_org[j];
    optflux[i][0] = optflux_org[j][0]/flux_mean_cont;
    optflux[i][1] = optflux_org[j][1]/flux_mean_cont;
  }
  ndrw = 1;
  date_span_cont = date_cont[nd_cont-1] - date_cont[0];
  

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
//        printf("%s\n", str);
        if(sscanf(str, "%lf %lf %lf\n", &date_line_org[ic], &hbb_org[ic][0], &hbb_org[ic][1]) < 3)
        {
          printf("# Wrong in reading %s.\n", parset.file_line);
          exit(0);
        }
//      printf("%lf %lf %lf\n", date[ic], optflux[ic][0], hbb[ic][0]);
        code_idx_line_org[ic] = idx;

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

    flux_mean_line = 0.0;
    idx = 0;
    for(i=0; i<obs_num_line[idx]; i++)
      flux_mean_line += hbb_org[i][0];
    flux_mean_line /= obs_num_line[idx];
    printf("Mean line flux of code %s:%f\n", code[idx], flux_mean_line);

  // sort over the data 
    gsl_sort_index(perm_line, date_line_org, 1, nd_line);
    for(i=0; i<nd_line; i++)
    {
      j = perm_line[i];
      date_line[i] = date_line_org[j];
      code_idx_line[i] = code_idx_line_org[j];
      hbb[i][0] = hbb_org[j][0]/flux_mean_line;
      hbb[i][1] = hbb_org[j][1]/flux_mean_line;
    }

    ndrw=2;
    date_span_line = date_line[nd_line-1] - date_cont[0];
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

    strcpy(tag[nt], "ScaleRangeLow");
    addr[nt] = &parset.scale_range_low;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ScaleRangeUp");
    addr[nt] = &parset.scale_range_up;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ShiftRangeLow");
    addr[nt] = &parset.shift_range_low;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ShiftRangeUp");
    addr[nt] = &parset.shift_range_up;
    id[nt++] = DOUBLE;

    
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
    parset.n_mcmc = 300000;
    parset.n_builtin = 100000;

    parset.scale_range_low = 0.5;
    parset.scale_range_up = 1.5;
    parset.shift_range_low = -1.0;
    parset.shift_range_up = 1.0;

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
  
  if(parset.n_mcmc < 50000 | parset.n_builtin < 50000)
  {
    printf("# n_mcmc and n_builtin should be larger than 50000. \n");
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
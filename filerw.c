#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

int memory_malloc()
{
  int i;

  list = malloc(nlist_max*sizeof(char *));
  for(i=0; i<nlist_max; i++)
  {
    list[i] = malloc(100*sizeof(char));
  }

// emission line  
  hbb = malloc(nlist_max * sizeof(double *));
  o3b = malloc(nlist_max * sizeof(double *));
  o3n = malloc(nlist_max * sizeof(double *));
  hbn = malloc(nlist_max * sizeof(double *));
  chi2 = malloc(nlist_max * sizeof(double));

  for(i=0; i<nlist_max; i++)
  {
    hbb[i] = malloc(10*sizeof(double));
    o3b[i] = malloc(6*sizeof(double));
    o3n[i] = malloc(6*sizeof(double));
    hbn[i] = malloc(4*sizeof(double));
  }

  o3flux = malloc(nlist_max * sizeof(double));
  o3center = malloc(nlist_max * sizeof(double));

// continuum
  optflux = malloc(nlist_max * sizeof(double *));
  optslope = malloc(nlist_max * sizeof(double *));
  chi2opt = malloc(nlist_max * sizeof(double));
  
  for(i=0; i<nlist_max; i++)
  {
    optflux[i] = malloc(2*sizeof(double));
    optslope[i] = malloc(2*sizeof(double));
  }


  date = malloc(nlist_max * sizeof(double));
  code = malloc(ncode_max * sizeof(char *));
  for(i=0; i<ncode_max; i++)
  {
    code[i] = malloc(100 * sizeof(char));
  }
  obs_num = malloc(ncode_max * sizeof(int));
  code_idx = malloc(nlist_max * sizeof(int));
  
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

int get_file_list(char *fname)
{
  FILE *fp;
  int i, eof;
  char str[100];

  fp=fopen(fname, "r");

  if(fp==NULL)
  {
    printf("Cannot open file %s\n", fname);
    return 1;
  }

  for(i=0; i<nlist_max; i++)
  {
    fgets(str, 100, fp);
    if(feof(fp))
      break;
    sscanf(str,"%s", list[i]);
//    printf("%d %s %d\n", i, list[i], feof(fp));
  }
  nlist = i;
  printf("There are %d files in total.\n", nlist);

  fclose(fp);
  return 0;
}

int get_o3flux()
{
  int i;

  for(i=0; i<nlist; i++)
  {
    o3flux[i] = o3b[i][0]+o3n[i][0];
  }

  return 0;
}

int get_o3center()
{
  const int nt=1000;
  double a1, a2, a3, ymax, xmax;
  double xl[nt], yl[nt], y[nt];
  int i, j;

  for(i=0; i<nt; i++)
  {
    xl[i] = 4900.0 + (5100.0-4900.0)/(nt-1) * i;
  }

  for(i=0; i<nlist; i++)
  {
    a1 = o3b[i][0];
    a2 = o3b[i][2];
    a3 = o3b[i][4];
    singlegauss(xl, yl, nt, a1, a2, a3); 

    for(j=0; j<nt; j++)
    {
      y[j] = yl[j];
    }

    a1 = o3n[i][0];
    a2 = o3n[i][2];
    a3 = o3n[i][4];
    singlegauss(xl, yl, nt, a1, a2, a3); 

    for(j=0; j<nt; j++)
    {
      y[j] += yl[j];
    }

    a1 = o3b[i][0] * flux1;
    a2 = o3b[i][2];
    a3 = o3b[i][4]-wave1;
    singlegauss(xl, yl, nt, a1, a2, a3); 

    for(j=0; j<nt; j++)
    {
      y[j] += yl[j];
    }

    a1 = o3n[i][0] * flux1;
    a2 = o3n[i][2];
    a3 = o3n[i][4]-wave1;
    singlegauss(xl, yl, nt, a1, a2, a3); 

    for(j=0; j<nt; j++)
    {
      y[j] += yl[j];
    }

    ymax = 0.0;
    for(j=0; j<nt; j++)
    {
      if(ymax<y[j])
      {
        xmax = xl[j];
        ymax = y[j];
      }
    }
    
    o3center[i] = xmax;
//    printf("%f\n", o3center[i]);
  }
}

int get_line_fitting_par(char *fname)
{
  FILE *flfp;
  int i, j, k, iter, eof;
  char str[100], fstr[300];

  flfp = fopen(fname, "r");
  if(flfp==NULL)
  {
    printf("cannot open file %s\n", fname);
  }

  for(j=0; j<nlist; j++)
  {
    fgets(fstr, 300, flfp);
    if(feof(flfp))
      break;
    sscanf(fstr, "%d %s %d %lf", &i, str, &iter, &chi2[j]);
    
    if(i!=j+1)
    {
      printf("i!=j in get_line\n");
      exit(-1);
    }
    
    if(strcmp(str, list[i-1]))
    {
      printf("wrong! in %s and %s\n", str, list[i-1]);
    }
//    printf("%d %s %s\n", i, str, list[i-1 ]);

// get broad Hbeta line
    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &hbb[i-1][0], &hbb[i-1][1], &hbb[i-1][2], &hbb[i-1][3]);
 //   printf("%f %f %f %f\n", hbb[i-1][0], hbb[i-1][1], hbb[i-1][2], hbb[i-1][3]);

    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &hbb[i-1][4], &hbb[i-1][5], &o3b[i-1][0], &o3b[i-1][1]);
//    printf("%f %f %f %f\n", hbb[i-1][4], hbb[i-1][5], o3b[i-1][0], o3b[i-1][1]);

    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &o3b[i-1][2], &o3b[i-1][3], &o3b[i-1][4], &o3b[i-1][5]);
//    printf("%f %f %f %f\n", o3b[i-1][2], o3b[i-1][3], o3b[i-1][4], o3b[i-1][5]);

    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &o3n[i-1][0], &o3n[i-1][1], &o3n[i-1][2], &o3n[i-1][3]);
//    printf("%f %f %f %f\n", o3n[i-1][0], o3n[i-1][1], o3n[i-1][2], o3n[i-1][3]);

    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &o3n[i-1][4], &o3n[i-1][5], &hbn[i-1][0], &hbn[i-1][1]);
//    printf("%f %f\n", o3n[i-1][4], o3n[i-1][5]);
    
    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf %lf %lf", &hbn[i-1][2], &hbn[i-1][3], &hbb[i-1][6], &hbb[i-1][7]);

    fgets(fstr, 300, flfp);
    sscanf(fstr, "%lf %lf", &hbn[i-1][8], &hbn[i-1][9]);

// read 11 times
    for(k=0; k<4; k++)
    {
      fgets(fstr, 300, flfp);
    }
  }

  fclose(flfp);
}


int get_cont_fitting_par(char *fname)
{
  FILE *fp;

  int i, j, k, iter, eof;
  char str[100], fstr[300];

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
    printf("cannot open file %s\n", fname);
  }
  
  for(j=0; j<nlist; j++)
  {
    fgets(fstr, 300, fp);
    if(feof(fp))
      break;
    sscanf(fstr, "%d %s %d %lf", &i, str, &iter, &chi2opt[j]);
    if(i!=j+1)
    {
      printf("i!=j in get_cont\n");
      exit(-1);
    }
//    printf("%d %s %d %lf", i, str, iter, chi2opt[i-1]);
    if(strcmp(str, list[i-1]))
    {
      printf("wrong! in %s and %s\n", str, list[i-1]);
    }
//    printf("%d %s %s\n", i, str, list[i-1 ]);


    fgets(fstr, 300, fp);

// get flux and slope    
    fgets(fstr, 300, fp);
    sscanf(fstr, "%lf %lf %lf %lf", &optflux[i-1][0], &optflux[i-1][1], &optslope[i-1][0], &optslope[i-1][1]);
//    printf("%f %f %f %f", optflux[i-1][0], optflux[i-1][1], optslope[i-1][0], optslope[i-1][1]);
// read 11 times
    for(k=0; k<2; k++)
    {
      fgets(fstr, 300, fp);
    }
  }
  fclose(fp);
}

void singlegauss(double *xin, double *y, int nd, double ain1, double ain2, double ain3)
{
  double a1, a2, a3, arg;
  double x;
  int i;

  a1 = ain1;
  a2 = ain2 * 1.0e5/C;
  a3 = log(ain3);

  for(i=0; i<nd; i++)
  {
    x = log(xin[i]);
    arg = (x - a3)/a2;
    y[i] = a1 *2.0*sqrt(log(2.0)/PI)/a2 * exp(-4.0*log(2.0)*arg*arg); 
  }

  for(i=0; i<nd; i++)
  {
    y[i] /= xin[i];
  }

}

void gausshermite(double *xin, double *y, int nd, double ain1, double ain2, 
  double ain3, double h3, double h4)
{
  double a1, a2, a3, arg, hx, f3, f4, h;
  double x;
  int i;


  a1 = ain1;
  a2 = ain2 * 1.0e5/C;
  a3 = log(ain3);

  for(i=0; i<nd; i++)
  {
    x = log(xin[i]);
    arg = (x - a3)/a2;
    hx=arg*(2.0*sqrt(2.0*log(2.0)));

    y[i] = a1 *2.0*sqrt(log(2.0)/PI)/a2 * exp(-4.0*log(2.0)*arg*arg); 
    f3 = 1.0/sqrt(6.0)*(2.0*sqrt(2.0)*hx*hx*hx-3.0*sqrt(2.0)*hx);
    f4=1.0/sqrt(24.0)*(4.0*hx*hx*hx*hx-12.0*hx*hx+3);
    h = 1.0 + h3*f3+h4*f4;
    y[i] *= h;
  }

  for(i=0; i<nd; i++)
  {
    y[i] /= xin[i];
  }

}

int scale_flux_by_OIII()
{
  double scale;
  int i;
  for(i=0; i<nlist; i++)
  {
    scale = o3flux_std / o3flux[i];
    optflux[i][0] *= scale;
    optflux[i][1] *= scale;
    hbb[i][0] *= scale;
    hbb[i][1] *= scale;
  }
  return 0;
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

  for(i=0; i<nlist; i++)
  {
    idx = code_idx[i];
    fprintf(fout, "%f\t%f\t%f\t%s\n", date[i], optflux[i][0] * ps_scale[idx] - es_scale[idx],
     optflux[i][1]* ps_scale[idx], code[code_idx[i]]);
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

  for(i=0; i<nlist; i++)
  {
    idx = code_idx[i];
    fprintf(fout, "%f\t%f\t%f\t%s\n", date[i], hbb[i][0] * ps_scale[idx], 
           hbb[i][1]* ps_scale[idx], code[code_idx[i]]);
  }

  fclose(fout);
}

int get_date_from_filename()
{
   int i, tmp;
   char str[100], *pstr;
	
   for(i=0; i<nlist; i++)
   {
     pstr = list[i];
     pstr = pstr + 2;
     strcpy(str, pstr);
     pstr = str + 5;
     *pstr = '\0';
     sscanf(str, "%d", &tmp);
//     printf("%d %s %d\n", i, str, tmp);
     date[i]=tmp;
//     printf("%s\t%f\n", str, date[i]);
   }
   return 0;
}

int read_dateset_code(char *fname)
{
  FILE *fin;
  int i, idx, ic;
  char str[100];
  
  fin=fopen(fname,"r");
  if(fin==NULL)
  {
    printf("cannot open file %s\n", fname);
    exit(-1);
  }
  
  for(i=0; i<ncode_max; i++)
  {
    fgets(str, 100, fin);
    if(feof(fin))
      break;
    sscanf(str, "%s %d", code[i], &obs_num[i]);
//    printf("%s %d\n", code[i], obs_num[i]);
  }
  ncode = i;
//  printf("%d\n", ncode);


// set the idx array
  idx = 0;
  ic = 0;
  for(i=0; i<nlist; i++)
  {
    ic++;
    code_idx[i] = idx;
//    printf("%d %d\n", i, code_idx[i]);
    if(ic==obs_num[idx])
    {
      ic=0;
      idx++;
    }
  }
  return 0;
}

int read_dataset(char *fname)
{
  FILE *fin;
  char str[300];
  int i, j, idx, ic;
  
  fin =fopen(fname, "r");
  if(fin==NULL)
  {
    printf("cannot open file %s\n", fname);
    exit(-1);
  }
  
  idx = 0;
  ic = 0;
  for(i=0; i<nlist_max; i++)
  {
    fgets(str, 300, fin);
    if(feof(fin))
      break;
    sscanf(str, "# %s %d\n", code[idx], &obs_num[idx]);
    printf("# %s %d\n", code[idx], obs_num[idx]);
    for(j=0; j<obs_num[idx]; j++)
    {
      fgets(str, 300, fin);
//      printf("%s\n", str);
      sscanf(str, "%lf %lf %lf %lf %lf", &date[ic], &optflux[ic][0], &optflux[ic][1], &hbb[ic][0], &hbb[ic][1]);
//      printf("%lf %lf %lf\n", date[ic], optflux[ic][0], hbb[ic][0]);
      code_idx[ic] = idx;
      
      optflux[ic][0] *= (o3flux_std);
      optflux[ic][1] *= (o3flux_std);
      hbb[ic][0] *= o3flux_std;
      hbb[ic][1] *= o3flux_std;
      
      //optflux[ic][1] = optflux[ic][0]*0.04;
      //hbb[ic][1] = hbb[ic][0]*0.035;
      
      //printf("%6.1f   %4.2f   %4.2f   %4.2f    %4.2f\n",  date[ic], optflux[ic][0], optflux[ic][1], hbb[ic][0], hbb[ic][1]);
      ic++;
    }
    idx++;
  }
  nlist = ic;
  ncode = idx;
  printf("%d %d\n", ic, idx); 
//  exit(-1);
}


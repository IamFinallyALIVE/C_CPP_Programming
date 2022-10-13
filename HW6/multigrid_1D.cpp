//author Peeyush Kumar 12th April,2021

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
  //#include <stdio.h>
  //#include <stdlib.h>
  //#include <math.h>
  //#include <complex.h>
using namespace std;

FILE * pFile = fopen ("solution.txt","w"); 
  
typedef struct{
    int N;
    int Lmax;
    int size[20];
    double a[20];
    double m;
    double scale[20];
  } param_t;

void relax(double *phi, double *res, int lev, int niter, param_t p);
void proj_res(double *res_c, double *rec_f, double *phi_f, int lev,param_t p);
void inter_add(double *phi_f, double *phi_c, int lev,param_t p);
double GetResRoot(double *phi, double *res, int lev, param_t p);

int main()
{  
  double *phi[20], *res[20];
  param_t p;
  int nlev;
  int i,lev;
  
  
  p.Lmax = 10;
  p.N = 2*(int)pow(2.,p.Lmax); 
  p.m = 0.001;
  nlev = 10; // NUMBER OF LEVELS
 
  if(nlev  > p.Lmax){ 
    printf("ERROR More levels than in lattice! \n");
    return 0; }
  
  printf("\n V cycle for  %d lattice nlev = %d out of max  %d \n \n", p.N, nlev, p.Lmax); 
  
  // initialize arrays__________________________________
  p.size[0] = p.N;
  p.a[0] = 1.0;
  p.scale[0] = 1.0/(2.0 + p.m*p.m);
  
  for(lev = 1;lev< p.Lmax+1; lev++) {
    p.size[lev] = p.size[lev-1]/2;
    p.a[lev] = 2.0 * p.a[lev-1];
    //  p.scale[lev] = 1.0/(2.0 + p.m*p.m*p.a[lev]*p.a[lev]);
    p.scale[lev] = 1.0/(2.0 + p.m*p.m);  // This seem to work better?? Why I am not sure!
  }


  for(lev = 0;lev< p.Lmax+1; lev++)
    {
      phi[lev] = (double *) malloc(p.size[lev] * sizeof(double)); //over ride with error to save space
      res[lev] = (double *) malloc(p.size[lev] * sizeof(double));
      for(i = 0;i< p.size[lev];i++)
       {
	  phi[lev][i] = 0.0;
          res[lev][i] = 0.0;
	};
    }  

   
  res[0][p.N/4] = 1.0*p.scale[0];  //Two unit point middle of fist half  N/2 and anti-source at reflected 
  res[0][3*p.N/4] = - 1.0*p.scale[0];
  
  // iterate to solve_____________________________________
  double resmag = 1.0; // not rescaled.
  int ncycle = 0; 
  int n_per_lev = 10;
  resmag = GetResRoot(phi[0],res[0],0,p);
  printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
 

  
  while(resmag > 0.000001 && ncycle < 100000 )
    { 
      ncycle +=1; 
      for(lev = 0;lev<nlev; lev++)   
	{    
	   relax(phi[lev],res[lev],lev, n_per_lev,p); 
	   proj_res(res[lev + 1], res[lev], phi[lev], lev,p);    
	}
      for(lev = nlev;lev >= 0; lev--) 
	{ 
  	  relax(phi[lev],res[lev],lev, n_per_lev,p);   
	  if(lev > 0) inter_add(phi[lev-1], phi[lev], lev, p);   
	}
      resmag = GetResRoot(phi[0],res[0],0,p);
      printf("At the %d cycle the mag residue is %g \n",ncycle,resmag);
    }


   int L  = p.size[0];
   for(int i = 0; i < p.N; i++)
     {
      double residual = res[0][i]/p.scale[0] - phi[0][i]/p.scale[0]  
	          + (phi[0][(i+1)%L] + phi[0][(i-1+L)%L] );
     fprintf(pFile, " %d     %g     %g      \n ",i,phi[0][i],residual);
     }
   
   fclose (pFile);
   
  return 0;
}

void relax(double *phi, double *res, int lev, int niter, param_t p)
{  
  int i, x,y;
   int L;
   L  = p.size[lev];
  
  for(i=0; i<niter; i++)    
    for(x = 0; x < L; x++)
	{
	  phi[x] = res[x] + p.scale[lev] * (phi[(x+1)%L] + phi[(x-1+L)%L] );			
       	}
  return;    
}

void proj_res(double *res_c, double *res_f, double *phi_f,int lev,param_t p)
{  
  int L, Lc, f_off, c_off, x, y;
  L = p.size[lev];
  double r[L]; // temp residue
  Lc = p.size[lev+1];  // course level
  
  //get residue
  for(x = 0; x< L; x++)
      r[x ] = res_f[x ] -  phi_f[x ] + p.scale[lev]*(phi_f[(x+1)%L] + phi_f[(x-1+L)%L]);
  
  //project residue
  for(x = 0; x< Lc; x++)  res_c[x ] = 0.5*( r[2*x]  + r[(2*x + 1)%L ] );

  return;
}

void inter_add(double *phi_f,double *phi_c,int lev,param_t p)
{ 
  int L, Lc, x, y;
  Lc = p.size[lev];  
  L = p.size[lev-1]; 
  for(x = 0; x< Lc; x++)
      {
	phi_f[2*x]              += phi_c[x];
	phi_f[(2*x + 1)%L ]     += phi_c[x];
      }
  //set to zero so phi = error 
  for(x = 0; x< Lc; x++) phi_c[x + y*L] = 0.0;
  return;
}

#if 0
double GetResRoot(double *phi, double *res, int lev, param_t p)
{ 
  int i, x;
  double residue;
  double ResRoot = 0.0;
  int L;
  if(lev != 0) cout<<"ERROR ON TOP LEVEL RES " << endl;
  L  = p.size[lev];
  for(x = 0; x < L; x++){
    residue = p.a[lev]*p.a[lev]*p.scale[lev]*b[x]-  phi[x]
            +  p.scale[0]*(phi[(x+1)%L] + phi[(x-1+L)%L]) ;
      ResRoot += residue*residue; 
    }

  return sqrt(ResRoot);    
}
#endif

double GetResRoot(double *phi, double *res, int lev, param_t p)
{ //true residue
  int i, x,y;
  double residue;
  double ResRoot = 0.0;
  int L;
  L  = p.size[lev];
  
  for(x = 0; x < L; x++) {
      residue = res[x]/p.scale[lev] - phi[x]/p.scale[lev]  
	+ (phi[(x+1)%L] + phi[(x-1+L)%L] );
      ResRoot += residue*residue;
    }

  return sqrt(ResRoot);    
}
#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <iomanip> 
#include <fstream>
using namespace std;
// 1D length
#define N 512

// Maximum number of iterations  100000000
#define ITER_MAX 50000

// How often to check the relative residual
#define RESID_FREQ 1000 

// The residual
#define RESID 1e-6

double magnitude(double T[N+2][N+2]);
void jacobi(double T[N+2][N+2], double b[N+2][N+2], double tmp[N+2][N+2]);
double getResid(double T[N+2][N+2], double b[N+2][N+2]);

int main(int argc, char** argv)
{
  int i,totiter,j;
  int done = 0;
  double T[N+2][N+2];
  double Ttmp[N+2][N+2]; 
  double b[N+2][N+2];
  double bmag, resmag;

  for (i=0;i<N+1;i++) {
      for (j=0;j<N+1;j++){ 
                T[i][j] = 0.0; 
                Ttmp[i][j] = 0.0; 
                b[i][j] = 0.0;
      }          
        }


//putting the hot source
  b[N/2][N/2] = 100.0;
  //Get magnitude of rhs
  bmag = magnitude(b);
  cout<<"bmag : "<<bmag<<endl;

  //std::ofstream myfile("jacobi_data.csv");
  for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
  {
    // do RESID_FREQ jacobi iterations
    jacobi(T, b, Ttmp);

    resmag = getResid(T, b);
    
    //myfile << totiter <<",";
    //myfile<< resmag/bmag <<","<<endl;
    cout<<totiter<<" "<<"residue : "<<resmag<<" bmag : "<<bmag<<" relative : "<<resmag/bmag<<endl;
    if (resmag/bmag < RESID) { done = 1; }
  }

 //myfile.close();

  return 0; 
}

double magnitude(double T[N+2][N+2])
{
   int i,j;
   double bmag;

   i = 0;
   j=0;
   bmag = 0.0;  
   for (i=0;i<N;i++) {
      for (j=0;j<N;j++){ 
            bmag = bmag + T[i][j]*T[i][j];
      }          
        }
   
   return sqrt(bmag);
}

void jacobi(double T[N+2][N+2], double  b[N+2][N+2], double  tmp[N+2][N+2])
{
  int iter,i,j;

  iter = 0; i = 0;j=0;

  for (iter=0;iter<RESID_FREQ;iter++)
  {
    #pragma omp for
    for (i=1;i<N;i++)
    {
       for (j=1;j<N;j++){ 

            tmp[i][j] = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j];
       }
    }


#pragma omp for
    for (i=1;i<N;i++)
    {
       for (j=1;j<N;j++){ 
              T[i][j] = tmp[i][j];
       }
    }

  }
}

double getResid(double  T[N+2][N+2], double b[N+2][N+2])
{
  int i,j;
  double localres,resmag;

  i = 0,j=0;
  localres = 0.0;
  resmag = 0.0;

  for (i=1;i<N;i++)
  {
    for (j=1;j<N;j++){ 
      localres = (b[i][j] - T[i][j] + 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) );
      localres = localres*localres;
      resmag = resmag + localres;
    }
  }

  resmag = sqrt(resmag);

  return resmag;
}


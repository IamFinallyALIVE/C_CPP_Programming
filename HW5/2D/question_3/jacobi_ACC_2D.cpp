#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace std;

// Maximum number of iterations
#define ITER_MAX 100000000

// How often to check the relative residual
#define RESID_FREQ 1000 

// The residual
#define RESID 1e-6

using namespace std;
int L,N;

double magnitude(double** T);
void jacobi(double** T, double** b, double** tmp);
double getResid(double** T, double** b);

int main(int argc, char** argv) {
  
  #pragma acc init
  for (N = 4; N<=10; N++) {
	L=pow(2,N);
  int i,j,totiter;
  double** T = new double*[L+1];
  double** Ttmp = new double*[L+1];
  double** b = new double*[L+1];
  double bmag, resmag;
  int done = 0;

  for (i=0;i<L+1;i++) { 
    T[i] = new double[L+1]; 
    Ttmp[i] = new double[L+1]; 
    b[i] = new double[L+1];  
    for (j = 0; j< L+1; j++) {
      T[i][j] = 0.0;
      Ttmp[i][j] = 0.0;
      b[i][j] = 0.0;
    }  }


  //update boundary condition
  for(i = 0; i<L + 1; i++){
    T[i][0] = sin(2* i * M_PI/(L + 1));
    T[i][L] = -sin(2* i* M_PI/(L + 1));
    b[i][0] = sin(2* i* M_PI/(L + 1));
    b[i][L] = -sin(2* i * M_PI/(L + 1));
  }
  
  for(j = 0; j< L + 1; j++){
    T[0][j] = sin(j * M_PI/(L + 1));
    T[L][j] = -sin(j * M_PI/(L + 1));
    b[0][j] = sin(j * M_PI/(L + 1));
    b[L][j] = -sin(j * M_PI/(L  +1));
  }
   
  //Get magnitude of rhs
  bmag = magnitude(b);

  // TIMING LINE 1: Get the starting timestamp. 
  std::chrono::time_point<std::chrono::steady_clock> begin_time = 
      std::chrono::steady_clock::now();

  for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
  {
    jacobi(T, b, Ttmp);

    resmag = getResid(T, b);

    if (resmag/bmag < RESID) { done = 1; }
  }

  // TIMING LINE 2: Get the ending timestamp.
  std::chrono::time_point<std::chrono::steady_clock> end_time =
      std::chrono::steady_clock::now();

   // TIMING LINE 3: Compute the difference.
  std::chrono::duration<double> difference_in_time = end_time - begin_time;

   // TIMING LINE 4: Get the difference in seconds.
  double difference_in_seconds = difference_in_time.count();

  printf("L: %d, %.10f seconds.\n", L, difference_in_seconds);

  for (int j =0; j<L+1;j++){
    delete[] T[j]; delete[] Ttmp[j]; delete[] b[j];
  }
  delete[] T; delete[] Ttmp; delete[] b;
}
  return 0;
}

double magnitude(double** T)
{
   int i;
   double bmag;
   bmag = 0.0;  
   for (i = 0; i<L+1; i++)
   {
    for (int j = 0; j<L+1; j++){
      bmag = bmag + T[i][j]*T[i][j];
   }}
   
   return sqrt(bmag);
}

void jacobi(double** T, double** b, double** tmp){
  
  for (int iter=0;iter<RESID_FREQ;iter++){
    #pragma acc parallel loop independent
    for (int i=1;i<L;i++){
      #pragma acc loop independent
      for (int j=1; j<L; j++)
        tmp[i][j] = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j];
    }

    #pragma acc parallel loop independent
    for (int i=1;i<L;i++){ 
      #pragma acc loop independent
      for (int j=1; j<L; j++)
        T[i][j] = tmp[i][j];
    }
  }
}

double getResid(double** T, double** b){
  int i;
  double localres,resmag;

  i = 0;
  localres = 0.0;
  resmag = 0.0;

  for (i=1;i<L;i++)
  {
    for (int j = 1; j<L; j++) {
      localres = 0.25*(b[i][j] + T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) - T[i][j];
      localres = localres*localres;
      resmag = resmag + localres;
    }
  }

  resmag = sqrt(resmag);

  return resmag;
}

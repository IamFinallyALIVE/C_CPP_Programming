//author Peeyush Kumar 12th April,2021
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <array>
#include <chrono>

using namespace std;
#include "fft.h"

int main()
{
int n;
for (n=1; n<=10;n++) {	
  int N = pow(2,n);
  Complex * omega = new Complex[N];
  Complex * F = new Complex[N];
  Complex * Fft1  = new Complex[N];
  Complex * Fnew = new Complex[N];
  Complex * Ftilde  = new Complex[N];
  
  makePhase(omega,N);

  /* Test slow FT */
  
  for(int x = 0; x < N; x++){
    F[x] = 2* sin( 2.0*PI*x/(double) N) + 4*cos( 2.0*PI*3.0*x/(double) N);

  }

  FFT1(F, Fft1, N); 
  
  clock_t FT_begin = clock();
  auto FT_begin_res = chrono::high_resolution_clock::now();
  FT(Ftilde, F,omega,N);
  auto FT_end_res = chrono::high_resolution_clock::now();
  auto FT_time_res = chrono::duration_cast<chrono::microseconds>(FT_end_res - FT_begin_res);
  clock_t FT_end = clock();
  double FT_time = 1000000.0 * (FT_end - FT_begin) / (double) CLOCKS_PER_SEC;

  
  clock_t FTinv_begin = clock();
  FTinv(Fnew,Ftilde, omega, N);
  clock_t FTinv_end = clock();
  double FTinv_time = 1000000.0 * (FTinv_end - FTinv_begin) / (double) CLOCKS_PER_SEC;
  
  
  clock_t FFT_begin = clock();
  auto FFT_begin_res = chrono::high_resolution_clock::now();
  FFT(Ftilde,Fft1, omega, N);
  auto FFT_end_res = chrono::high_resolution_clock::now();
  auto FFT_time_res = chrono::duration_cast<chrono::microseconds>(FFT_end_res - FFT_begin_res);
  clock_t FFT_end = clock();
  double FFT_time = 1000000.0 * (FFT_end - FFT_begin) / (double) CLOCKS_PER_SEC;
  
  std::cout << std::fixed <<std::setprecision(20) << N <<" "<< FT_time <<" "<<  FFT_time << std::endl;
}

  return 0;
}
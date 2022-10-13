#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <iostream>
using namespace std;
//#define N 512

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 1000

// The residual
#define RESID 1e-2

// Useful globals
int world_size; // number of processes
int my_rank; // my process number
int L,N;


double magnitude(double** T, const int size);
void jacobi(double** T, double** b, double** Ttmp, const int size);
double getResid(double** T, double** b, const int size);

int main(int argc, char** argv)
{
   int i,j,totiter;
   int done = 0;
   //double **T, **Ttmp, **b; 
   double bmag, resmag;
   int local_size;
   
   // Initialize MPI
   MPI_Init(&argc, &argv);
   
   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   
   // Get the rank
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   

   // N is L in this case
   for (N = 4; N<=10; N++) {
   done = 0;
	L=pow(2,N);
   // lx and ly is one fourth
   int Lx = L/4;
   int Ly = L/4;

   double** T = new double* [Lx];
   double** Ttmp = new double* [Lx];
   double** b = new double* [Lx];

   for (i=0;i<Lx;i++) {
      T[i] = new double[Ly]; 
      Ttmp[i] = new double[Ly]; 
      b[i] = new double[Ly]; 
      for(int j = 0; j < Ly; j++){
         T[i][j] = 0.0;
         Ttmp[i][j] = 0.0;
         b[i][j] = 0.0; 
      }
  }
// Special Boundary conditions as mentioned in HW5
   if (my_rank%Ly == 0) {
      int Lxx = my_rank/Lx;
      for(i = 0; i<Lx; i++){
         T[0][i] = sin(2*(Lxx*Lx+i)*M_PI/(Lx*Ly+1));
         b[0][i] = sin(2*(Lxx*Lx+i)*M_PI/(Lx*Ly+1));
      }
   }

   if (my_rank%Ly == Ly-1) {
      int Lxx = my_rank/Lx;
      for(i = 0; i<Lx; i++){
         T[Ly-1][i] = -sin(2*(Lxx*Lx+i)*M_PI/(Lx*Ly+1));
         b[Ly-1][i] = -sin(2*(Lxx*Lx+i)*M_PI/(Lx*Ly+1));
      }
   }

   if (my_rank/Ly == 0) {
      int Lyy = my_rank;
      for(j = 0; j<Ly; j++){
         T[j][0] = sin((Lyy*Ly+j)*M_PI/(Ly*Lx+1));     
         b[j][0] = sin((Lyy*Ly+j)*M_PI/(Ly*Lx+1));
     }
   }

   if (my_rank/Ly == Lx-1){
      int Lyy = my_rank % Ly;
      for(j = 0; j<Ly; j++){
       T[j][Lx-1] = -sin((Lyy*Ly+j)*M_PI/(Ly*Lx+1));
       b[j][Lx-1] = -sin((Lyy*Ly+j)*M_PI/(Ly*Lx+1));
      }
   }

   bmag = magnitude(b, Lx);


   // TIMING LINE 1: Get the starting timestamp. 
   std::chrono::time_point<std::chrono::steady_clock> begin_time = 
      std::chrono::steady_clock::now();
   
   for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
   {
  
      jacobi(T, b, Ttmp, Lx);
      resmag = getResid(T, b, Lx);
      
      if (resmag/bmag < RESID) { done = 1; }
   }

   std::chrono::time_point<std::chrono::steady_clock> end_time =
      std::chrono::steady_clock::now();
   std::chrono::duration<double> difference_in_time = end_time - begin_time;
   double seconds = difference_in_time.count();
   
   if (my_rank == 0) {
      printf("%d %.10f \n", L, seconds);
   }
   

   for (i = 0; i < Lx; i++) {
      delete[] T[i]; delete[] Ttmp[i]; delete [] b[i];
   }

   delete[] T; delete[] Ttmp; delete [] b;
   }
  
   MPI_Finalize();
   
   return 0;
}


double magnitude(double** T, const int size)
{
   int i;
   double bmag = 0.0;
   double global_bmag = 0.0 ; // used for global reduce!

   for (i = 0; i < size; i++)
   {
      for(int j = 0; j < size; j++)
         bmag = bmag + T[i][j]*T[i][j];
   }
   
   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}


void jacobi(double** T, double** b, double** Ttmp, const int size){
   int iter,i, j, k;
   
   // Prepare for async send/recv
   MPI_Request request[8];
   int requests;
   MPI_Status status[8];
   
   const int lower_limit = (my_rank % size == 0 ) ? 1 : 0;
   const int upper_limit = (my_rank % size == size-1) ? size-1 : size;
   const int top_limit= (my_rank/size == 0) ? 1 : 0;
   const int bot_limit = (my_rank/size == size-1) ? size-1 : size;

   iter = 0; i = 0;

   double *left_buffer = new double[size];
   double *right_buffer = new double[size];

   double *top_buffer = new double[size];
   double *bottom_buffer = new double[size];

   for (int j = 0; j< size; j++) {
      left_buffer[j] = 0.0;
      right_buffer[j] = 0.0;
   }
   for (int j = 0; j< size; j++) {
      top_buffer[j] = 0.0;
      bottom_buffer[j] = 0.0;
   }
   
   for (iter=0;iter<RESID_FREQ;iter++) {

         requests=0;
      
         // Fill the left buffer. Send to the right, listen from the left.
         MPI_Isend(&T[size-1][0],  size, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&left_buffer[0], size, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         
         // Fill the right buffer. Send to the left, listen from the right.
         MPI_Isend(&T[0][0], size, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&right_buffer[0], size, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);

         
         double* top_Ttmp = new double[size];
         double* bot_Ttmp = new double[size];
         for (k=0; k<size; k++) {
            top_Ttmp[k] = T[k][0]; bot_Ttmp[k] = T[k][size-1];
         }

         MPI_Isend(&top_Ttmp[0], size, MPI_DOUBLE, (my_rank+world_size-N/16)%world_size, 2, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&bottom_buffer[0], size, MPI_DOUBLE, (my_rank+N/16)%world_size, 2, MPI_COMM_WORLD, request + requests++);
		 
		 
         MPI_Isend(&bot_Ttmp[0],  size, MPI_DOUBLE, (my_rank+N/16)%world_size, 3, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&top_buffer[0], size, MPI_DOUBLE, (my_rank+world_size-N/16)%world_size, 3, MPI_COMM_WORLD, request + requests++);
        
         // Loop over the rest.
         for (i=1;i<size-1;i++)
         {
            for (int j = 1; j<size-1; j++){
               Ttmp[i][j] = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j];
            } 
         }

         // Wait for async.
         MPI_Waitall ( requests, request, status );
         
         // BCs
         if (my_rank % size != 0)
         {
            i = 0;
            for (int j = 1; j < size-1; j++)
               Ttmp[i][j] = 0.25*(T[i+1][j]+left_buffer[j]+T[i][j+1]+T[i][j-1]) + b[i][j];
         }

         if (my_rank % size != size-1)
         {
            i = size -1;
            for (int j = 1; j < size-1; j++)
               Ttmp[i][j] = 0.25*(right_buffer[j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j];
         }

         if (my_rank / size != 0)
         {
            j = 0;
            for (i = 1; i < size-1; i++)
               Ttmp[i][j] = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+top_buffer[i]) + b[i][j];
         }

         if (my_rank / size != size-1)
         {
            j =  size -1;
            for (i = 1; i < size-1; i++)
               Ttmp[i][j] = 0.25*(T[i+1][j]+T[i-1][j]+bottom_buffer[i]+T[i][j-1]) + b[i][j];
         }

         // update BCs
         if(lower_limit == 0 && top_limit== 0)
            Ttmp[0][0] = 0.25*(T[1][0]+left_buffer[0]+T[0][1]+top_buffer[0]) + b[0][0];
         if(upper_limit == size && top_limit== 0)
            Ttmp[size-1][0] = 0.25*(right_buffer[0]+T[size-2][0]+T[size-1][1]+top_buffer[size-1]) + b[size-1][0];
         if(lower_limit == 0 && bot_limit == size)
            Ttmp[0][size-1] = 0.25*(T[1][size-1]+left_buffer[size-1]+bottom_buffer[0]+T[0][size-2]) + b[0][size-1];
         if(upper_limit == size && bot_limit == size)
            Ttmp[size-1][size-1] = 0.25*(right_buffer[size-1]+T[size-2][size-1]+bottom_buffer[size-1]+T[size-1][size-2]) + b[size-1][size-1];
         
         for (i =lower_limit;i<upper_limit;i++)
         {
            for (j = top_limit; j < bot_limit; j++)
               T[i][j] = Ttmp[i][j];
         }
      }

   delete[] left_buffer; delete[] right_buffer; delete[] top_buffer; delete[] bottom_buffer;
   
   MPI_Barrier(MPI_COMM_WORLD);
}

double getResid(double** T, double** b, const int size){
   int i;
   double localres,resmag;
   double global_resmag;
  
   MPI_Request request[8];
   int requests;
   MPI_Status status[8];

   const int lower_limit = (my_rank % size == 0 ) ? 1 : 0;
   const int upper_limit = (my_rank % size == size-1) ? size-1 : size;
   const int top_limit= (my_rank/size == 0) ? 1 : 0;
   const int bot_limit = (my_rank/size == size-1) ? size-1 : size;
   
   // grab the left and right buffer.
   double *left_buffer = new double[size];
   double *right_buffer = new double[size];

   // grab the top and bottom buffer.
   double *top_buffer = new double[size];
   double *bottom_buffer = new double[size];

   // init buffer
   for (int j = 0; j< size; j++) {
      left_buffer[j] = 0.0;
      right_buffer[j] = 0.0;
   }
   for (int j = 0; j< size; j++) {
      top_buffer[j] = 0.0;
      bottom_buffer[j] = 0.0;
   }
   
   requests=0;


   // Fill the left buffer. Send to the right, listen from the left.
   MPI_Isend(&T[size-1][0],  size, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&left_buffer[0], size, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

   // Fill the right buffer. Send to the left, listen from the right.
   MPI_Isend(&T[0][0], size, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&right_buffer[0], size, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);

   
   double* top_Ttmp = new double[size];
   double* bot_Ttmp = new double[size];
   for (int k=0; k<size; k++) {
      top_Ttmp[k] = T[k][0];
      bot_Ttmp[k] = T[k][size-1];
   }

   // Fill the bottom buffer. Send to the top, listen from top
   MPI_Isend(&bot_Ttmp[0],  size, MPI_DOUBLE, (my_rank+4)%world_size, 3, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&top_buffer[0], size, MPI_DOUBLE, (my_rank+world_size-4)%world_size, 3, MPI_COMM_WORLD, request + requests++);

   // Fill the top buffer. Send to the bottom, listen from the bottom.
   MPI_Isend(&top_Ttmp[0], size, MPI_DOUBLE, (my_rank+world_size-4)%world_size, 2, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&bottom_buffer[0], size, MPI_DOUBLE, (my_rank+4)%world_size, 2, MPI_COMM_WORLD, request + requests++);
   
   i = 0;
   localres = 0.0;
   global_resmag = 0.0;
   resmag = 0.0;

   // Loop over rest.
   for (i=1;i<size-1;i++)
   {
      for (int j = 1; j < size-1; j++) {
         localres = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j] - T[i][j];
         localres = localres*localres;
         resmag = resmag + localres;
      } 
   }
 
   // Wait for async.
   MPI_Waitall ( requests, request, status );
   
   if (my_rank % size != 0)
   {
      i = 0;  
      for (int j = 1; j < size-1; j++){
         localres = 0.25*(T[i+1][j]+left_buffer[j]+T[i][j+1]+T[i][j-1]) + b[i][j] - T[i][j];
         localres = localres*localres;
         resmag = resmag + localres;
      }

   }

   if (my_rank % size != size-1)
   {
      i = size -1;
      for (int j = 1 ; j < size-1; j++){
         localres = 0.25*(right_buffer[j]+T[i-1][j]+T[i][j+1]+T[i][j-1]) + b[i][j] - T[i][j];
         localres = localres*localres;
         resmag = resmag + localres;
      }
   }

   if (my_rank / size != 0)
   {
      int j = 0;
      for (i = 1; i < size-1; i++){
          localres = 0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+top_buffer[i]) + b[i][j] - T[i][j];
         localres = localres*localres;
         resmag = resmag + localres;
      }
   }

   if (my_rank / size != size-1)
   {
      int j =  size -1;
      for (i = 1; i < size-1; i++){
         localres = 0.25*(T[i+1][j]+T[i-1][j]+bottom_buffer[i]+T[i][j-1]) + b[i][j] - T[i][j];
         localres = localres*localres;
         resmag = resmag + localres;
      }
   }
	// update BCs
   if(lower_limit == 0 && top_limit== 0){
      localres = 0.25*(T[1][0]+left_buffer[0]+T[0][1]+top_buffer[0]) + b[0][0]- T[0][0];
      localres = localres*localres;
      resmag = resmag + localres;
   }
   if(upper_limit == size && top_limit== 0){
      localres = 0.25*(right_buffer[0]+T[size-2][0]+T[size-1][1]+top_buffer[size-1]) + b[size-1][0]- T[size-1][0];
      localres = localres*localres;
      resmag = resmag + localres;
   }
   if(lower_limit == 0 && bot_limit == size){
      localres = 0.25*(T[1][size-1]+left_buffer[size-1]+bottom_buffer[0]+T[0][size-2]) + b[0][size-1]- T[0][size-1];
      localres = localres*localres;
      resmag = resmag + localres;
   }

   if(upper_limit == size && bot_limit == size){
      localres = 0.25*(right_buffer[size-1]+T[size-2][size-1]+bottom_buffer[size-1]+T[size-1][size-2]) + b[size-1][size-1]- T[size-1][size-1];
      localres = localres*localres;
      resmag = resmag + localres;
   }

   delete[] left_buffer; delete[] right_buffer; delete[] top_buffer; delete[] bottom_buffer;
   // Reduce. 
   MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

   return sqrt(global_resmag);
}

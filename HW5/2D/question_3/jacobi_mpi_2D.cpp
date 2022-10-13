#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <chrono>

using namespace std;

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 10000

// The residual
#define RESID 1e-2

// Useful globals
int world_size; // number of processes
int my_rank; // my process number
int N, L;

double magnitude(double** T, const int size, const int L);
void jacobi(double** T, double** b, double** Ttmp, const int size, const int L);
double getResid(double** T, double** b, const int size, const int L);

int main(int argc, char** argv)
{
	int totiter;
	int done = 0;
	double bmag, resmag;
	int local_size; 
   
   // Initialize MPI
   MPI_Init(&argc, &argv);

   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
      
   // Get the rank
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   for(int N = 4; N <= 10; N++)
   {
	  L=pow(2,N);
      local_size = L/world_size;
      if (my_rank == (world_size-1)) { local_size += (L % world_size) ; }

	  int i,j,totiter;
      double **T = new double *[local_size];
      double **Ttmp = new double *[local_size];
      double **b = new double *[local_size];
      for (i = 0; i < local_size; i++)
      {
            T[i] = new double[L + 1];
            Ttmp[i] = new double[L + 1];
            b[i] = new double[L + 1];
      }

      for (i = 0; i < local_size; i++)
      {
            for (j = 0; j < L + 1; j++)
            {
                  T[i][j] = 0.0;
                  Ttmp[i][j] = 0.0;
                  b[i][j] = 0.0;
            }
      }

      int source_rank = (L/2)/(L/world_size);
      if (my_rank == source_rank) { b[L / 2 - source_rank * local_size][L / 2] = 1.0; }

      bmag = magnitude(b, local_size, L);
       
	  // TIMING LINE 1: Get the starting timestamp. 
	  std::chrono::time_point<std::chrono::steady_clock> begin_time = std::chrono::steady_clock::now();
	  
      for (totiter = RESID_FREQ; totiter < ITER_MAX && done == 0; totiter += RESID_FREQ)
      {
         jacobi(T, b, Ttmp, local_size, L);        

         resmag = getResid(T, b, local_size, L);
         
         if (resmag/bmag < RESID) 
         {
            done = 1; 
         } }
		std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();
		std::chrono::duration<double> difference_in_time = end_time - begin_time;
		double seconds = difference_in_time.count();
			
      if (my_rank == 0) {
         printf("%d %.10f\n", L, seconds);
      }
	  
      for (i = 0; i < local_size; i++)
      {
            delete[] T[i], Ttmp[i], b[i];
      }
      delete[] T; delete[] Ttmp; delete[] b;
   }

   MPI_Finalize();
   
   return 0;
}

double magnitude(double** T, const int size, const int L)
{
   double bmag = 0.0;
   double global_bmag = 0.0;
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size - 1) ? size - 1 : size;

   for (int i = lower_limit; i < upper_limit; i++)
   {
      for (int j = 1; j < L; j++)
      {
            bmag = bmag + T[i][j] * T[i][j];
      }
   }
   
   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}

void jacobi(double** T, double** b, double** Ttmp, const int size, const int L)
{  
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;
   
   // grab the left and right buffer.
   double *left_buffer = new double[L + 1];
   double *right_buffer = new double[L + 1];

   {
      for (int iter = 0; iter < RESID_FREQ; iter++)
      {

         requests = 0;
      
         // Fill the left buffer. Send to the right, listen from the left.
         MPI_Isend(T[size-1], L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(left_buffer, L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

         // Fill the right buffer. Send to the left, listen from the right.
         MPI_Isend(T[0],   L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(right_buffer, L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);
         
         for (int i = 1; i < size - 1; i++)
         {
               for (int j = 1; j < L; j++)
               {
                     Ttmp[i][j] = 0.25 * (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) + b[i][j];
               }
         }

         // Wait for async.
         MPI_Waitall ( requests, request, status );
         
         
         if (my_rank != 0)
         {
            for (int i = 1; i < L; i++)
            {
               Ttmp[0][i] = 0.25 * (T[1][i] + left_buffer[i] + T[0][i + 1] + T[0][i - 1]) + b[0][i];
            }
         }

         if (my_rank != world_size - 1)
         {
            for (int i = 1; i < L; i++)
            {
               Ttmp[size - 1][i] = 0.25 * (T[size - 2][i] + right_buffer[i] + T[size - 1][i + 1] + T[size - 1][i - 1]) + b[size - 1][i];
            }     }

         for (int i = lower_limit; i < upper_limit; i++)
         {
            for (int j = 1; j < L; j++)
            {
                  T[i][j] = Ttmp[i][j];
            }
         } } }
   
   MPI_Barrier(MPI_COMM_WORLD);

   delete [] left_buffer;
   delete [] right_buffer;
}

double getResid(double** T, double** b, const int size, const int L)
{
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size - 1) ? size - 1 : size;

   double localres,resmag;
   double global_resmag;
   
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   // grab the left and right buffer.
   double *left_buffer = new double[L + 1];
   double *right_buffer = new double[L + 1];
   
   requests=0;
      
   // Fill the left buffer. Send to the right, listen from the left.
   MPI_Isend(T[size-1], L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(left_buffer, L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

   // Fill the right buffer. Send to the left, listen from the right.
   MPI_Isend(T[0],   L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(right_buffer, L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);

   localres = 0.0;
   global_resmag = 0.0;
   resmag = 0.0;

   // Loop over rest.
   for (int i = 1; i < size - 1; i++)
   {
      for (int j = 1; j < L; j++)
      {
            localres = T[i][j] - (0.25 * (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) + b[i][j]);
            localres = localres * localres;
            resmag = resmag + localres;
      }
   }
   
   // Wait for async.
   MPI_Waitall ( requests, request, status );
   
   // Impose zero bc.
   if (my_rank != 0)
   {
      for (int i = 1; i < L; i++)
      {
         localres = T[0][i] - (0.25 * (T[1][i] + left_buffer[i] + T[0][i + 1] + T[0][i - 1]) + b[0][i]);
         localres = localres * localres;
         resmag = resmag + localres;
      }
   }

   // Impose zero bc.
   if (my_rank != world_size - 1)
   {
      for (int i = 1; i < L; i++)
      {
         localres = T[size - 1][i] - (0.25 * (T[size - 2][i] + right_buffer[i] + T[size - 1][i + 1] + T[size - 1][i - 1]) + b[size - 1][i]);
      }
   }


   // Reduce. 
   MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

   delete [] left_buffer; delete [] right_buffer;
   
   return sqrt(global_resmag);
}

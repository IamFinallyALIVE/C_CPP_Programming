#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <iostream>

using namespace std;
//#define N 512 --- needs to go from 64 to 16384

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 1000

// The residual
#define RESID 1e-2 

// Useful globals
int world_size; // number of processes
int my_rank; // my process number



double magnitude(double** x, int L, const int size);
void jacobi(double** x, double** b, int L, double** tmp, const int size);
double getResid(double** x, double** b, int L, const int size);

int main(int argc, char** argv)
{

   FILE *outfile;
   outfile = fopen("jacobi_mpi_2d.dat","w");


   
   // Initialize MPI
   MPI_Init(&argc, &argv);
   
   // Get the number of processes
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   
   // Get the rank
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   
   // this loop willl vary L as mentoned in HW5
   for (int L = 16; L<=512; L=L*2) {
   
   int i,j,totiter;
   int done = 0;
   double **Temp, **Ttmp, **b;    
   double bmag, resmag, mag;
   int local_size;

      // Figure out my local size. The last rank gets the leftover. 
      local_size = L/world_size;
      
      if (my_rank == (world_size-1)) { local_size += (L % world_size) ; }
      
      Temp = new double*[local_size];
      Ttmp = new double*[local_size];
      b = new double*[local_size];
	  
	for (i=0;i<local_size;i++) { 
	 Temp[i] = new double[L + 1];
	 Ttmp[i] = new double[L + 1];
	 b[i] = new double[L + 1];
	}

      for (i=0;i<local_size;i++) { 
         for(j=0; j<L+1; j++) {Temp[i][j] = 0.0; Ttmp[i][j] = 0.0; b[i][j] = 0.0; }}
      
      // The source only lives on a particular rank!
	  int source_rank = (L/2)/(L/world_size);

      if (my_rank == source_rank) { b[L/2 - source_rank*(L/world_size)][L/2 - source_rank*(L/world_size)] = 100.0; }
      
      //Get magnitude of rhs
      bmag = magnitude(b, L, local_size);
      //printf("bmag: %.8e\n", bmag);

      // TIMING LINE 1: Get the starting timestamp.
      std::chrono::time_point<std::chrono::steady_clock> begin_time = std::chrono::steady_clock::now();

      for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
      {

         // do RESID_FREQ jacobi iterations
         jacobi(Temp, b, L, Ttmp, local_size);

         getResid(Temp, b, L, local_size);
                
         if (resmag/bmag < RESID) { done = 1; }
      }


      std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();

      std::chrono::duration<double> difference_in_time = end_time - begin_time;

      double seconds = difference_in_time.count();

	if (my_rank==0) {
      fprintf(outfile,"%d %.8e\n", L, seconds);
      cout << "L " << L << " Time_In_Seconds " << seconds << endl;
	}
      delete[] Temp; delete[] Ttmp; delete[] b;
}
   
   fclose(outfile);
   MPI_Finalize();
   
   return 0;
}

double magnitude(double** x, int L, const int size)
{
   int i,j;
   double bmag;
   double global_bmag; // used for global reduce!
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;

   i = 0; j=0;
   bmag = 0.0;  
   global_bmag = 0.0;
   for (i = lower_limit; i<upper_limit; i++)
   {
	   for (j = 1; j<L; j++)
     bmag = bmag + x[i][j]*x[i][j];
   }
   

   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}

void jacobi(double** Temp, double** b, int L, double** tmp, const int size)
{
   int iter,i,j,k;
   double value;
   
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;
   
   // grab the left and right buffer.
   double* left_buffer = new double[L + 1];
   double* right_buffer = new double[L + 1];
   

   //iter = 0; i = 0;

      for (iter=0;iter<RESID_FREQ;iter++)
      {
         double mag = 0.;

         requests=0;
      
         // Fill the left buffer. Send to the right, listen from the left.
         MPI_Isend(Temp[size-1],   L+1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(left_buffer, L+1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

         // Fill the right buffer. Send to the left, listen from the right.
         MPI_Isend(Temp[0],   L+1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(right_buffer, L+1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
        
         // Loop over the rest.
         for(i=1; i<size-1;i++ ) {
         for (j=1;j<L;j++)
         {
           tmp[i][j] = b[i][j] + 0.25 * (Temp[i+1][j]+Temp[i-1][j]+Temp[i][j+1]+Temp[i][j-1]);
         }
		 }
         // Wait for async.
         MPI_Waitall ( requests, request, status );
         
         
         // Impose zero bc on one side
         if (my_rank != 0)
         {
			for(i=1; i<L;i++ ) {
            tmp[0][i] = 0.25 * (Temp[0][i+1]+Temp[0][i-1]+Temp[1][i]+left_buffer[i]) + b[0][i];
         }
		 }
         
         // Impose zero bc on other side
         if (my_rank != world_size-1)
         {
			for(i=1; i<L;i++ ) { 
            tmp[size-1][i] = 0.25 * (Temp[size-1][i+1]+Temp[size-1][i-1]+Temp[size-2][i]+right_buffer[i]) + b[size-1][i];
         }
         } 
		for(i=lower_limit; i<upper_limit;i++ ) {
         for (k=1;k<L;k++)
         {
            Temp[i][k] = tmp[i][k];
         }
		}
   }
   
   MPI_Barrier(MPI_COMM_WORLD);
}


double getResid(double** Temp, double** b,  int L, const int size)
{
   int i,j;
   double localres,resmag;
   double global_resmag;
  
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   // grab the left and right buffer.
   double* left_buffer = new double[L + 1];
   double* right_buffer = new double[L + 1];
   
   requests=0;
      
   // Fill the left buffer. Send to the right, listen from the left.
   MPI_Isend(Temp[size-1],   L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(left_buffer, L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

   // Fill the right buffer. Send to the left, listen from the right.
   MPI_Isend(Temp[0], L + 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(right_buffer,L + 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);

   i = 0;
   localres = 0.0;
   global_resmag = 0.0;
   resmag = 0.0;

   // Loop over rest.
   for (i=1;i<size-1;i++)
   {
	 for (j=1;j<L;j++){   
      localres = Temp[i][j] - (b[i][j] + 0.25 * (Temp[i+1][j]+Temp[i-1][j]+Temp[i][j+1]+Temp[i][j-1]));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   }
   
   // Wait for async.
   MPI_Waitall ( requests, request, status );
   
   // Impose zero bc.
   if (my_rank != 0)
   {
	for (i=1;i<L;i++){ 
      localres = Temp[0][i]-(0.25 * (Temp[0][i+1]+Temp[0][i-1]+Temp[1][i]+left_buffer[i]) + b[0][i]);
      localres = localres*localres;
      resmag = resmag + localres;
   }
   }

   if (my_rank != world_size-1)
   {
	for (i=1;i<L;i++){
      localres = Temp[size-1][i]-(0.25 * (Temp[size-1][i+1]+Temp[size-1][i-1]+Temp[size-2][i]+right_buffer[i]) + b[size-1][i]);
      localres = localres*localres;
      resmag = resmag + localres;
   }
   }
   // Reduce. 
   MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_resmag);
}

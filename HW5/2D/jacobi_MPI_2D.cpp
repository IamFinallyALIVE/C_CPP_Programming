#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <iostream>

using namespace std;

// Maximum number of iterations
#define ITER_MAX 10000000

// How often to check the relative residual
#define RESID_FREQ 1000

// The residual
#define RESID 1e-6 // as in early problem

// Useful globals
int world_size; // number of processes
int my_rank; // my process number



double magnitude(double* x, const int size);
double jacobi(double* x, double* b, int L, double* tmp, const int size);
double getResid(double* x, double* b, const int size);

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
   
   for (int L = 16; L<512; L=L*2) {
   
   int i,totiter;
   int done = 0;
   double **Temp, **Ttmp, **b;     // replace previous xtmp etc make them two-dimensional
   double bmag, resmag, mag;
   int local_size;

      // Figure out my local size. The last rank gets the leftover. 
      local_size = L/world_size;
      
      if (my_rank == (world_size-1)) { local_size += (L % world_size) ; }

      //printf("I am rank %d of %d and I have a local size %d.\n", my_rank, world_size, local_size); 
      
      Temp = new double[L+2];
      Ttmp = new double[L+2];
      b = new double[L+2];

      for (i=0;i<L+2;i++) { 
         for(j=0; j<L+2; j++) {Temp[i][j] = 0.0; Ttmp[i][j] = 0.0; b[i][j] = 0.0; }}
      
      // b[N/2] = 1.0;
      // The source only lives on a particular rank!
      b = [L/2][L/2];
      // TIMING LINE 1: Get the starting timestamp.
      std::chrono::time_point<std::chrono::steady_clock> begin_time = std::chrono::steady_clock::now();

      for (totiter=RESID_FREQ;totiter<ITER_MAX && done==0;totiter+=RESID_FREQ)
      {

         // do RESID_FREQ jacobi iterations
         mag = jacobi(Temp, b, L, Ttmp, local_size);

         //resmag = getResid(x, b, local_size);
                
         if (mag < RESID) { done = 1; }
      }

      std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();
      std::chrono::duration<double> difference_in_time = end_time - begin_time;
      double seconds = difference_in_time.count();

      fprintf(outfile,"%d %.8e\n", L, seconds);
      cout << "L " << L << " difference_in_seconds " << seconds << endl;

      free(Temp); free(Ttmp); free(b);
}
   
   fclose(outfile);
   // Clean up
   MPI_Finalize();
   
   return 0;
}

double magnitude(double* x, const int size)
{
   int i;
   double bmag;
   double global_bmag; // used for global reduce!
   const int lower_limit = (my_rank == 0) ? 1 : 0;
   const int upper_limit = (my_rank == world_size-1) ? size-1 : size;

   i = 0;
   bmag = 0.0;  
   global_bmag = 0.0;
   for (i = lower_limit; i<upper_limit; i++)
   {
     bmag = bmag + x[i]*x[i];
   }


   // Reduce. 
   MPI_Allreduce(&bmag, &global_bmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_bmag);
}

double jacobi(double* Temp, double* b, int L, double* tmp, const int size)
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
   double left_buffer = 0.0;
   double right_buffer = 0.0;
   

   //iter = 0; i = 0;

      for (iter=0;iter<RESID_FREQ;iter++)
      {
         double mag = 0.;
         for(i=1; i<L;i++ ) {
         requests=0;
      
         // Fill the left buffer. Send to the right, listen from the left.
         MPI_Isend(&Temp[size-1],   1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

         //printf("I am rank %d of %d and I received %.8e from the left.\n", my_rank, world_size, left_buffer);
         
         // Fill the right buffer. Send to the left, listen from the right.
         MPI_Isend(&Temp[0],   1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
         //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);
         
        
         // Loop over the rest.
         for (j=1;j<size-1;j++)
         {
           tmp[i][j] = b[i][j] + 0.25 * (Temp[i+1][j]+Temp[i-1][j]+Temp[i][j+1]+Temp[i][j-1]);
           mag += (Temp[i][j]-tmp[i][j]) * (Temp[i][j]-tmp[i][j]);
         }

         // Wait for async.
         MPI_Waitall ( requests, request, status );
         
         
         // Impose zero bc on one side
         if (my_rank != 0)
         {
            tmp[i][0] = 0.25 * (Temp[i+1][0]+Temp[i-1][0]+Temp[i][1]+left_buffer) + b[x][0];
         }
         
         // Impose zero bc on other side
         if (my_rank != world_size-1)
         {
            tmp[i][size-1] = 0.25 * (Temp[i+1][size-1]+Temp[i-1][size-1]+Temp[i][size-2]+right_buffer) + b[x][size-1];
            0.5*(right_buffer+x[size-2]) + b[size-1];
         }
          

         for (k=lower_limit;k<upper_limit;k++)
         {
            Temp[i][k] = tmp[i][k];
         }
      }
   value = mag;
   }
   return value;
   //MPI_Barrier(MPI_COMM_WORLD);
}

double getResid(double* x, double* b, const int size)
{
   int i;
   double localres,resmag;
   double global_resmag;
  
   // Prepare for async send/recv
   MPI_Request request[4];
   int requests;
   MPI_Status status[4];
   
   // grab the left and right buffer.
   double left_buffer = 0.0;
   double right_buffer = 0.0;
   
   requests=0;
      
   // Fill the left buffer. Send to the right, listen from the left.
   MPI_Isend(&x[size-1],   1, MPI_DOUBLE, (my_rank+1)%world_size, 1, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&left_buffer, 1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 1, MPI_COMM_WORLD, request + requests++);

   //printf("I am rank %d of %d and I received %.8e from the left.\n", my_rank, world_size, left_buffer);

   // Fill the right buffer. Send to the left, listen from the right.
   MPI_Isend(&x[0],   1, MPI_DOUBLE, (my_rank+world_size-1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   MPI_Irecv(&right_buffer, 1, MPI_DOUBLE, (my_rank+1)%world_size, 0, MPI_COMM_WORLD, request + requests++);
   //printf("I am rank %d of %d and I received %.8e from the right.\n", my_rank, world_size, right_buffer);

   i = 0;
   localres = 0.0;
   global_resmag = 0.0;
   resmag = 0.0;

   // Loop over rest.
   for (i=1;i<size-1;i++)
   {
      localres = (b[i] - x[i] + 0.5*(x[i+1] + x[i-1]));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   
   // Wait for async.
   MPI_Waitall ( requests, request, status );
   
   // Handle boundaries, acknowledging 0 bcs.
   // Impose zero bc.
   if (my_rank != 0)
   {
      localres = (b[0] - x[0] + 0.5*(x[1] + left_buffer));
      localres = localres*localres;
      resmag = resmag + localres;
   }
   // Impose zero bc.
   if (my_rank != world_size-1)
   {
      localres = (b[size-1] - x[size-1] + 0.5*(right_buffer + x[size-2]));
      localres = localres*localres;
      resmag = resmag + localres;
   }

   // Reduce. 
   MPI_Allreduce(&resmag, &global_resmag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
   
   return sqrt(global_resmag);
}

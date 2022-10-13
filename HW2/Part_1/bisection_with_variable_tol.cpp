#include <iostream>
#include <iomanip>
#include <cmath>
#include <float.h> // For long doubles.
#include <fstream>

using std::cout;
using std::cin;
using std::endl;
using std::ios;

using std::abs;

// This function performs one iteration of Newton's method
// and returns a new guess (x - f(x)/f'(x) -> x_new).
// For now, you need to hard-code the numerical function and
// its derivative.
long double newton(long double x, long double A, double N);

// This function performs one iteration of bisection
// and updates either "min" or "max" (note how they are both
// passed by reference), and returns the current "midpoint".
// Again, you need to hard-code the numerical function. Bisection
// does not require a derivative.
long double bisection(long double A, double N, long double & min, long double & max);



int main()
{
    
    std::ofstream myfile("bisection_data.csv");
    double tol ;
    for(int N=2;N<4;N++){
        cout<<"N = "<<N<<endl;
    for(int i =1; i<16;i++){
        // Declare a tolerance
        tol = pow(10.0,-i);
        
        
        
        // Declare variables to hold the current guess
        // and relative error.
        long double x = 0.0, fractional_error = 0.5;
        long A=2;
        // Declare a counter.
        int count;
        
        // Fix the output format.
        cout.setf(ios::fixed,ios::floatfield);
        cout.precision(20);
        cout<< " Bisection Method starting with min = 0 and max = A\n";
        long double min, max;
        min = 0.0;
        max = A;
        count = 0;

        do
        {
          count++;
          x =  bisection( A, N,  min, max);
          fractional_error = 0.5*abs(pow(x,N)/A-1);
          //cout.precision(20);
          //cout << x   << "\t" <<  fractional_error << endl;
        }
        while(fractional_error > tol);
        cout.precision(40);
        
        
        cout<<"iteration : "<<count<<" tol :"<<tol<<" ans : "<<x<<" error : "<<(x - pow(A, 1.0/N))<<endl;
        
        myfile << count <<",";
        myfile<< tol <<","<<endl;
    }
    }
    myfile.close();
  return  0;
}


// This routine is currently hard coded for the function
// f(x) - x^2 - A
long double bisection(long double A, double N, long double & min, long double & max)
{
  long double x  = (min + max)/2.0;
  if(pow(x,N)-A < 0.0)
    min = x;
  else
    max = x;

  return x;
}

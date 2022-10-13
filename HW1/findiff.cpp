#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

# define DATA_SAVE_FLAG 1

using namespace std;

double f(double x)
{
  return sin(x);
}

double derivative(double x)
{
  return cos(x);
}

double forward_diff(double x, double h)
{
  return (f(x+h)-f(x))/h;
}

double backward_diff(double x, double h)
{
  return (f(x)-f(x-h))/h;
}

double central_diff(double x, double h)
{
  return (f(x+0.5*h)-f(x-0.5*h))/h;
}

double relative_error(double diff, double exact){
    return (diff-exact)/(exact);
};


int main(int argc, char** argv)
{
  int i;
  double h;
  const double x = 1.0;

  cout << setprecision(15); // .. 10 points after the decimal.
    
    //csv for h
    // This line initiates the csv file creation
   

    
        std::ofstream myfile("data.csv");
  
  // Loop over 17 values of h: 1 to 10^(-16).
  h = 1.0;
    cout<<"h  "<<"                  Sin(x)  "<<"                 Cos(x)  "<<"        Forward Diff  "<<"      Central Diff : "<<"     Backward : "<<endl;
    
  
  for (i = 0; i < 17; i++)
  {
      h=h - 0.055;
      
      cout<<h<<"    "<<sin(x)<<"    "<<derivative(x)<<"    "<<forward_diff(x, h)<<"    "<<central_diff(x, h)<<"    "<<backward_diff(x, h)<<endl;
      
      //ordr of values in csv : h, backward_error, forward error, central_error,
      
 
      myfile << h <<",";
      myfile<< relative_error(backward_diff(x, h), derivative(x)) <<",";
      myfile<< relative_error(forward_diff(x, h), derivative(x))<<"," ;
      myfile<< relative_error(central_diff(x, h), derivative(x))<<endl;
      
      
      
    // Print stuff to graph
  }
   
        myfile.close();
    
  return 0;
}




    std::ofstream myfile_x_8("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/x_8_data.csv");
    std::ofstream myfile_cos("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/cos_data.csv");
    std::ofstream myfile_func_2("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/func_2_data.csv");

    double x_8_exact=0.2222222222222222;
    double cosine_exact=1.273239544735163;
    double func_2_exact=1.570796326794897;

     for (order = 2; order <= 32; order=order+2){
      
      myfile_x_8 << order <<",";
      myfile_x_8<< error(x_8_exact, gauss_quad_integral(order,x_8)) <<",";
      myfile_x_8<< error(x_8_exact, Trapzoidal_rule(order,x_8)) <<",";

      myfile_cos << order <<",";
      myfile_cos<< error(cosine_exact, gauss_quad_integral(order,cosine)) <<",";
      myfile_cos<< error(cosine_exact, Trapzoidal_rule(order,cosine)) <<",";

      myfile_func_2 << order <<",";
      myfile_func_2<< error(func_2_exact, gauss_quad_integral(order,func_2)) <<",";
      myfile_func_2<< error(func_2_exact, Trapzoidal_rule(order,func_2)) <<",";

      
      
      
    // Print stuff to graph
  }
   
        myfile.close();
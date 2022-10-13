//
//  main.cpp
//  practice
//
//  Created by Peeyush Kumar on 2/25/22.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <ios>
#include <fstream>
#define PI 3.14159265
#define   SIZE   32
using namespace std;

double gauss_elimination( double a[SIZE][SIZE],double b[SIZE],int n,std::vector<double> v,double (*func)(double)){

    float x[SIZE], ratio;
    int i,j,k;
    cout<< setprecision(15)<< fixed;

    for(i=1;i<=n;i++)
     {
        a[i][n+1]=b[i];
     }
    

    // perform Gauss Elimination
     for(i=1;i<=n-1;i++)
     {
          if(a[i][i] == 0.0)
          {
               cout<<"Failed";
               exit(0);
          }
          for(j=i+1;j<=n;j++)
          {
               ratio = a[j][i]/a[i][i];

               for(k=1;k<=n+1;k++)
               {
                      a[j][k] = a[j][k] - ratio*a[i][k];
               }
          }
     }
     // extract values
     x[n] = a[n][n+1]/a[n][n];

     for(i=n-1;i>=1;i--)
     {
          x[i] = a[i][n+1];
          for(j=i+1;j<=n;j++)
          {
                  x[i] = x[i] - a[i][j]*x[j];
          }
          x[i] = x[i]/a[i][i];
     }

     // print solution
     cout<< endl<<"Solution: "<< endl;
     double ans=0;
     for(i=1;i<=n;i++)
     {
          cout<<"x["<< i<<"] = "<< v[i-1]<<"  "<<"w["<< i<<"] = "<< x[i]<< endl;
          ans=ans+ x[i]*func(v[i-1]);
     }

    std::cout<<"Integral using guass : "<<ans<<std::endl;
    return ans;

}

std::vector<double> getLegendreCoeff(int m)
{
    if (m == 0)
    {
        return {1};
    }
    if (m == 1)
    {
        return {0, 1};
    }

    std::vector<double> coeffs(m + 1);

    
    std::vector<double> v =  getLegendreCoeff(m - 1);
    std::vector<double> u =  getLegendreCoeff(m - 2);

    
    double a = (2.0 * m - 1.0) / m;
    double b = (m - 1.0) / m;

    int first = 1;
    if ( m % 2 == 0 )
    {
        coeffs[0] = -b * u[0];
        first = 2;
    }
    for (int i = first; i < m - 1; i += 2)
    {
        coeffs[i] = (a * v[i - 1] - b * u[i]);
    }
    coeffs[m] = a * v[m - 1];

    return coeffs;
}

float legendre_poly(float n,float x){
    if(n==0){
        return float(1.0);
    }
    if (n==1){ return float(x);}

    return float((((2*n)-1)/n)*x*legendre_poly(n-1,x) - ((((n)-1)/n)*legendre_poly(n-2,x)));
    
}

// forgot to multiply derivative with constants
double compute_function(std::vector<double> coeff,double x){
    long length=coeff.size();
    double sum=0;
    for (int i=0;i<=length;i++){
        sum= sum+coeff[i]*pow(x, i);
    }
    return sum;
}

double compute_function_derivative(std::vector<double> coeff,double x){
    long length=coeff.size();
    double sum=0;
    for (int i=0;i<=length;i++){
        sum= sum+ (i+1)*coeff[i]*pow(x, i);
    }
    return sum;
}

void print_vector(std::vector<double>v){
    long length=v.size();
    for (int i = 0;i<length;i++){
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
}
std::vector<double> find_roots(std::vector<double>poly,std::vector<double>derivative,double n,double tolerance, int tolerance_for_dis){
    std::vector<double> roots;
    // initiliaze roots with formula given with equation 9 of HW2
    long length=poly.size();
    
    for (int i = 1;i<length;i++){
        roots.push_back(double((1 - (1/(8*n*n))+ (1/(8*n*n*n)))*cos(PI*((4*i -1)/(4*n +2)))));
    }
    
    double x_0;
    double x_new, displacement ;
    int dis_tolerance=0;
    double old_displacement =0;
    for (int i = 0;i<length-1;i++){
        x_0=roots[i];
        do{
            x_new=x_0-(compute_function(poly, x_0)/(compute_function_derivative(derivative, x_0) ));
            displacement=abs(abs(x_new)-abs(x_0));
            x_0=x_new;
            
            // to kep chek if the displacement is varying
            
            if (abs(old_displacement-displacement) == 0){
                dis_tolerance=dis_tolerance+1;
            } else{
                dis_tolerance=0;
            }
            old_displacement=displacement;
            if (dis_tolerance>=tolerance_for_dis){
                break;
            }
            //break;
        }while(displacement>tolerance);
        roots[i]=x_new;
        //std::cout<<std::endl;
    }
    
    //std::cout<<"The roots are as follows : "<<std::endl;
    
    return roots;
}

std::vector<double>legend_driver(int n){
    if(n>33){
        std::cout<<"ERROR : n is greater than 33"<<'\n';
    }
    std::vector<double> coef =getLegendreCoeff(n);
    
    std::vector<double> derivative;
        long length=coef.size();
    for (int i=1;i<length;i++){
        derivative.push_back(coef[i]);
    }
    
    return find_roots(coef, derivative, n, 0.000001,10);
    
    
}


double x_8(double x){

    return pow(x,8);
}

double cosine(double x){

        return  cos(x*(PI/2));
}

double func_2(double x){
        return (1/(x*x +1));
}

double Trapzoidal_rule(double n, double (*func)(double) ){

        //double h= ((1-(-1))/n);
        double a,b;  
        a=-1;
        b=1;
        double h= (b-a)/n;
        double f_x0=func(a);
        double f_xn=func(b);

        double ans= f_x0 + f_xn;
        for(int i=1;i<n;i++){
            ans=ans+ 2*(func(a+i*h));
        }   

        std::cout<<"Integral using Trap : "<<ans*(h/2)<<std::endl;

        return ans*(h/2);

}

double gauss_quad_integral( double order, double (*func)(double)  ){
    int i,j;
    std::vector<double> a=legend_driver(order);
    double A[SIZE][SIZE], b[SIZE];
   
 
    //initialize the A matrix with legendre polynomals
    
    for(i=1;i<=order;i++)
         {
              for(j=1;j<=order;j++)
              {
                   A[i][j]= pow(a[j-1], i-1) ;
              }
         }
    
    //initialize the b matrix with zeros;
    //load vector into array b
    for(i=1;i<=order;i++)
         {
             b[i]= double((1- pow((-1),i) )/i) ;
         }
    
     return gauss_elimination(A,b,order,a,func);
}

double error(double exact_value, double calulated_value){

    return abs((calulated_value-exact_value)/exact_value);
   
}
int main()
{
    
    //legend_driver(); 
    int order,i,j;
    std::cout << "What order Legendre polynomial? ";
    std::cin >> order;
    gauss_quad_integral(order,func_2);
    Trapzoidal_rule(order,func_2);
    
    /*
     std::ofstream myfile_x_8("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/x_8_data.csv");
    std::ofstream myfile_cos("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/cos_data.csv");
    std::ofstream myfile_func_2("/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/func_2_data.csv");

    double x_8_exact=0.2222222222222222;
    double cosine_exact=1.273239544735163;
    double func_2_exact=1.570796326794897;

     for (order = 2; order <= 30; order=order+2){
      
      myfile_x_8 << order <<",";
      myfile_x_8<< error(x_8_exact, gauss_quad_integral(order,x_8)) <<",";
      myfile_x_8<< error(x_8_exact, Trapzoidal_rule(order,x_8)) <<std::endl;

      myfile_cos << order <<",";
      myfile_cos<< error(cosine_exact, gauss_quad_integral(order,cosine)) <<",";
      myfile_cos<< error(cosine_exact, Trapzoidal_rule(order,cosine)) <<std::endl;

      myfile_func_2 << order <<",";
      myfile_func_2<< error(func_2_exact, gauss_quad_integral(order,func_2)) <<",";
      myfile_func_2<< error(func_2_exact, Trapzoidal_rule(order,func_2)) <<std::endl;

  }
   
        myfile_x_8.close();
        myfile_cos.close();
        myfile_func_2.close();
*/
    return 0;
}








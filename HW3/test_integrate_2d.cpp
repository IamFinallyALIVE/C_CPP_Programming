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
#define   SIZE   34
using namespace std;

double gauss_elimination( double a[SIZE][SIZE],double b[SIZE],int n,std::vector<double> v,double (*func)(double,double)){

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
     double ans=0;
     for(i=1;i<=n;i++)
     {

          //cout<<"x["<< i<<"] = "<< v[i-1]<<"  "<<"w["<< i<<"] = "<< x[i]<< endl;
          for(j=1;j<=n;j++){
              ans=ans+ x[j]*x[i]*func(v[i-1],v[j-1]);
          }
     }

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


double exp_cos_sin(double x,double y){

    return exp(-pow(x,2) - (pow(y,2))/8) * cos(PI*x) *sin((PI/8)*y);
    // exact value 0
}

double sqrt_x_y(double x,double y){

        return  sqrt(x*x - y*y + 2);
        // exact value 5.62421
}

double func_x8_y8(double x,double y){

        return ( pow(x,8) + pow(y,8) +(pow((y-1),3)*pow((x-3),5)));
        // exact value = 2688.88888889
}


double gauss_quad_integral( double order, double (*func)(double,double)  ){
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
    int order;
    //std::cout << "What order Legendre polynomial? ";

    for(order =2;order<=32;order=order+2){

            std::cout<<"Gauss Integral for N = "<<order<<" for sqrt_x_y is  "<<gauss_quad_integral(order,sqrt_x_y)<<std::endl;
            
            std::cout<<"Gauss Integral for N = "<<order<<" for func_x8_y8 is  "<<gauss_quad_integral(order,func_x8_y8)<<std::endl;
            
            std::cout<<"Gauss Integral for N = "<<order<<" for exp_cos_sin is  "<<gauss_quad_integral(order,exp_cos_sin)<<std::endl;

            std::cout<<endl;
    }
    return 0;
}








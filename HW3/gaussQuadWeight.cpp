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
#define PI 3.14159265
#define   SIZE   32
using namespace std;

void gauss_elimination( double a[SIZE][SIZE],double b[SIZE],int n,std::vector<double> v){

    float x[SIZE], ratio;
    int i,j,k;
    cout<< setprecision(10)<< fixed;

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
     for(i=1;i<=n;i++)
     {
          cout<<"x["<< i<<"] = "<< v[i-1]<<"  "<<"w["<< i<<"] = "<< x[i]<< endl;
     }
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

void print_coeff(std::vector<double> const& v);

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
        std::cout<<std::endl;
    }
    
    //std::cout<<"The roots are as follows : "<<std::endl;
    
    return roots;
}

std::vector<double>legend_driver(int n){
    if(n>30){
        std::cout<<"ERROR : n is greater than 30"<<'\n';
    }
    std::vector<double> coef =getLegendreCoeff(n);
    
    std::vector<double> derivative;
        long length=coef.size();
    for (int i=1;i<length;i++){
        derivative.push_back(coef[i]);
    }
    
    return find_roots(coef, derivative, n, 0.000001,10);
    
    
}


int main()
{
    
    //legend_driver();
    int order,i,j;
    std::cout << "What order Legendre polynomial? ";
    std::cin >> order;
    
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
    
    gauss_elimination(A,b,order,a);
    return 0;
}


void print_coeff(std::vector<double> const& v)
{
    long i=v.size();
    while (i>=0){
        i=i-1;
        if(i>0){
            std::cout<< v[i]<<"x^"<<i<<" + ";
        }
        else{
            std::cout<< v[0]<<"x^"<<i;
            break;
        }
        
    }
    std::cout<<'\n';
}







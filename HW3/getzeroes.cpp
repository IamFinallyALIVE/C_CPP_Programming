//
//  main.cpp
//  practice
//
//  Created by Peeyush Kumar on 2/17/22.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <ios>
#define PI 3.14159265

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
    //print_vector(roots);
    
    return roots;
}

int main()
{
    int n;
    std::cout<<"What order Legendre polynomial?"<<std::endl;
    std::cin>>n;
    if(n>30){
        std::cout<<"ERROR : n is greater than 30"<<'\n';
        return 0;
    }
    std::vector<double> coef =getLegendreCoeff(n);
    
    std::vector<double> derivative;
        long length=coef.size();
    for (int i=1;i<length;i++){
        derivative.push_back(coef[i]);
    }
    
    find_roots(coef, derivative, n, 0.000001,10);
    
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






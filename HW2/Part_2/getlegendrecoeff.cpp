//
//  main.cpp
//  practice
//
//  Created by Peeyush Kumar on 2/17/22.
//

#include <iostream>
#include <vector>
#include <cmath>


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

    // Initialize with zero, only at most (half + 1) of the terms will be changed later
    std::vector<double> coeffs(m + 1);

    // Consider some form of memoization instead of this recursion
    std::vector<double> v =  getLegendreCoeff(m - 1);
    std::vector<double> u =  getLegendreCoeff(m - 2);

    // using literals of floating point type, 'm' is promoted by the compiler
    double a = (2.0 * m - 1.0) / m;
    double b = (m - 1.0) / m;

    int first = 1;
    // If 'm' is even, so is (m - 2) and its first element is zero. It can be skipped.
    // It also avoids annoying '-0' in the output
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


int main()
{
    int n;
    std::cout<<"What order Legendre polynomial?"<<std::endl;
    std::cin>>n;
    print_coeff(getLegendreCoeff(n));
    std::cout<<'\n';
    
    
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

}






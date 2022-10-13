   
#include<iostream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>

#define   SIZE   32

using namespace std;


void gauss_elimination( double a[SIZE][SIZE],double b[SIZE],int n){

    float x[SIZE], ratio;
	int i,j,k;
    cout<< setprecision(3)<< fixed;

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
	  	cout<<"x["<< i<<"] = "<< x[i]<< endl;
	 }
}

int main()
{
	 double a[SIZE][SIZE],b[SIZE];
	 int i,j,n;

    cout<< setprecision(3)<< fixed;
	 cout<<"Enter number of unknowns: ";
	 cin>>n;
	 for(i=1;i<=n;i++)
	 {
		  for(j=1;j<=n;j++)
		  {
			   a[i][j]= double(pow((-1 +(2*(j-1)/(n-1))),(i-1))) ;
               cout<<a[i][j]<<endl;
		  }
	 }
	 cout<<"Enter labels "<< endl;

     for(j=1;j<=n;j++){
         b[j]= double((1-pow(-1,i))/i) ;
     }


     cout<<"original answer_2"<<endl;
     gauss_elimination(a,b,n);
	
	 return(0);
}

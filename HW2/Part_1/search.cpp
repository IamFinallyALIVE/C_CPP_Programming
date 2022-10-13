#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "search.h"
#include "sort.h"

using std::cout;
using std::endl;

// A function to print arrays -- debug sorting for small arrays
void printArray(int list[], int arraySize)
{
  cout << " --------" << endl;
  for (int i = 0; i < arraySize; i++)
  {
    cout << list[i] <<  " " << endl;
  }
  cout << " --------" << endl;
}


void perform_diff_searches(int Nsize, double p[6]){
    
    int * list = new  int[Nsize];

    // Iterators.
    int i;

    // Counters
    int BinCount; // counter for binary search.
    int LinCount; // counter for linear search.
    int DctCount; // counter for dictionary search.
       
    int time_seconds = time(0);
    srand(time_seconds%100); // seed the random number generator.

    // Fill the array with random numbers.
    for(i = 0; i < Nsize; i++)
    {
      list[i] =  rand()%Nsize;
    }

    // Randomly choose an element to search for.
    int find = list[rand()%Nsize] ;
    i = linearSearch(list,find ,Nsize, LinCount);
    //Perform a binary search.
    i = binarySearch(list,find ,Nsize, BinCount);
    // Perform a dictionary search.
    i = dictionarySearch(list,find ,Nsize, DctCount);
    // Free the allocated memory.
    delete [] list;
    p[0]=double(Nsize/2);
    p[1]=log(Nsize)/log(2);
    p[2]= log(log(Nsize)/log(2))/log(2);
    p[3]=float(LinCount);
    p[4]=double(BinCount);
    p[5]=double(DctCount);
   
}

int main()
{
 
    double p[6];
    std::ofstream myfile("search_timing.dat");
    std::ofstream my_csv_file("data.csv");
    
    int N=0;
    for (int i =1;i<8;i++){
        cout<<"Value of i "<<i<<endl;
        N=pow(10, i);
        perform_diff_searches(N,p);
        
        cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<p[3]<<" "<<p[4]<<" "<<p[5]<<endl;
        
        myfile << N <<",";
        myfile<< p[3] <<",";
        myfile<< p[4]<<"," ;
        myfile<< p[5]<<endl;
        
        my_csv_file << N <<",";
        my_csv_file<< p[0] <<",";
        my_csv_file<< p[1]<<"," ;
        my_csv_file<< p[2]<<endl;
        
    }
    myfile.close();
  return 0;
}



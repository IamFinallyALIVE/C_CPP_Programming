#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 19:01:09 2022

@author: peeyushkumar
"""



import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

''' in tolerance and iterations 0-14 index is for N=2 and 15-end is for N=3'''
iterations=[]
tolerance=[]



#['h,Backwad', 'Relative', 'Error,Forward', 'Relative', 'Error,Central', 'Relatie', 'Error']

''' Please Enter Directory here for the data '''
directory='/Users/peeyushkumar/Desktop/Parallel Prog/HW2_code/bisection_data.csv'

with open(directory, newline='') as csvfile:
    
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for line in spamreader:
        temp=line[0].split(',')

        iterations.append(float(temp[0]))
        tolerance.append(float(temp[1]))



def plot_data(x,y,N, curve_flag=False,Curve=None,popt_n_2=None,popt_n_3=None):
    plt.scatter(x,y)
    #plt.scatter(x2,y2)
    plt.xlabel('Tolerance Value ---->')
    plt.ylabel('Iterations Count ---->')
    plt.title('Tolerance vs Iterations '+N)
    
    if curve_flag!=False:
        if N=='Square root':
            plt.plot(Curve(np.array(x),popt_n_2[0],popt_n_2[1],popt_n_2[2]),y)
        elif N=='Cube root':
            plt.plot(Curve(np.array(x),popt_n_3[0],popt_n_3[1],popt_n_3[2]),y)
    plt.show()
    

def plot_complexities():
    
    N=np.array([10**x for x in range(8)])
   
    plt.plot(N,N/2.0)
    plt.plot(N,np.log(N))
    plt.plot(N,np.log(np.log(N)))
    plt.xlabel('N ---->')
    plt.ylim([0,90])
    plt.ylabel('Complexities Value---->')
    plt.title('N vs Complexities')
    plt.legend(['N/2','log(N)','log(log(N))'])
    plt.show()
    
def curve(N,c0,c1,c2):
    return c0 + (c1*np.log(1/N)) +(c2*np.log(np.log(1/N)))
   
    
''' fitting curve for square root'''
popt_n_2, _ = curve_fit(curve, np.array(tolerance[:15]),np.array(iterations[:15]) )

''' fitting curve for cube root'''
popt_n_3, _ = curve_fit(curve,  np.array(tolerance[15:]),np.array(iterations[15:]))



''' Plot graph for N=2'''
plot_data(tolerance[:15],iterations[:15],'Square root')
''' Plot graph for N=3'''
plot_data(tolerance[15:],iterations[15:],'Cube root')  


'''plot data and fitted curve' for N=2'''
plot_data(tolerance[:15],iterations[:15],'Square root',curve_flag=True,Curve=curve,popt_n_2=popt_n_2)


'''plot data and fitted curve for N=3'''
plot_data(tolerance[15:],iterations[15:],'Cube root',curve_flag=True,Curve=curve,popt_n_3=popt_n_3)

'''plotiing complexities'''
plot_complexities()



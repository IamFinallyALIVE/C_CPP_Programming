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



#['h,Backwad', 'Relative', 'Error,Forward', 'Relative', 'Error,Central', 'Relatie', 'Error']

''' Please Enter Directory here for the data '''
files=[{'file':"cos_data.csv", 'name':'Cosine((pi*x)/2)' },
       {'file':"func_2_data.csv", 'name':'1/(x^2 + 1)' },
       {'file':"x_8_data.csv", 'name':'x^8' }]

directory='/Users/peeyushkumar/Desktop/Parallel_Prog/HW3code/'


def plot_data(x,y1,y2,title=None,xlabel='N--->',ylabel='Relative Error'):
    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.legend(['Gauss Relative Error','Trap Relative Error'])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Graph showing Relative Error for '+(title))
    plt.show()
    

for i in range(len(files)):
    with open(directory+files[i]['file'], newline='') as csvfile:
        n=[]
        gauss=[]
        trap=[]
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for line in spamreader:
            temp=line[0].split(',')
    
            
            n.append(float(temp[0]))
            gauss.append(float(temp[1]))
            trap.append(float(temp[2]))
        
        plot_data(n,gauss,trap,title=files[i]['name'])




    

''' Plot graph for N=2'''
#plot_data(tolerance[:15],iterations[:15],'Square root')
''' Plot graph for N=3'''
#plot_data(tolerance[15:],iterations[15:],'Cube root')  



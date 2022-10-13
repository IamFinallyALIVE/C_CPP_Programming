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
data=[]
res=[]



#['h,Backwad', 'Relative', 'Error,Forward', 'Relative', 'Error,Central', 'Relatie', 'Error']

''' Please Enter Directory here for the data '''
directory=['gauss_seidel_data','jacobi_data','red_black_data']

for file in directory:
    with open((file+'.csv'), newline='') as csvfile:
        temp_data=[]
        temp_res=[]
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for line in spamreader:
            temp=line[0].split(',')
    
            temp_data.append(float(temp[0]))
            temp_res.append(float(temp[1]))
        data.append(temp_data)
        res.append(temp_res)



def plot_data(x,y):
    
    for X,Y in zip(x,y):
        plt.plot(X,Y)
        
    plt.xlabel('Iterations---->')
    plt.ylabel('Relative Residue---->')
    plt.legend(['gauss_seidel','jacobi','red_black'])
    plt.show()
    
plot_data(data,res)

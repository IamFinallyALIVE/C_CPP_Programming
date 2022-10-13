#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 09:12:07 2022

@author: peeyushkumar
"""



import csv
import matplotlib.pyplot as plt


h=[]
backward_error=[]
forward_error=[]
central_error=[]



#['h,Backwad', 'Relative', 'Error,Forward', 'Relative', 'Error,Central', 'Relatie', 'Error']

''' Please Enter Directory here for the data '''
directory='data.csv'

with open(directory, newline='') as csvfile:
    
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for line in spamreader:
        temp=line[0].split(',')
        h.append(float(temp[0]))
        backward_error.append(float(temp[1]))
        forward_error.append(float(temp[2]))
        central_error.append(float(temp[3]))



'''plotting the graph'''

plt.plot(h,backward_error)
plt.plot(h,forward_error)
plt.plot(h,central_error)
plt.title('h vs Error Trade-off')
plt.xlabel('h---->')
plt.ylabel('Error----->')
plt.legend(['Backward_error','Forward_error','Central_error'])
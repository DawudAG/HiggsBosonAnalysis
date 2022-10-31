#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 19:45:25 2021

@author: dawudabd-alghani
"""

import numpy as np
import matplotlib.pyplot as plt
import STOM_higgs_tools as sht
import scipy as sp
import scipy.optimize

def loadData():
    f = plt.figure(1)
    vals = sht.generate_data()
    bin_heights, bin_edges, patches = plt.hist(vals, range = [104.0,155.0], bins = 150)
    centerVal = []
    for i in range(0,len(bin_edges)-1):
        midVal = (bin_edges[i]+bin_edges[i+1])*0.5
        centerVal.append(midVal)
    return vals,centerVal,bin_heights,bin_edges

def one():
    vals,centerVal,bin_heights,bin_edges = loadData()

    for i in range(0, len(bin_heights)-1):
        print(bin_heights[i])
        print(bin_edges[i+1] - bin_edges[i])
    std = np.sqrt(bin_heights)
    plt.errorbar(centerVal, bin_heights, yerr=std, fmt='r.')
    plt.show()
    input()


def twoa():
    total = 0
    length = 0 
    vals,centerVal,bin_heights,bin_edges = loadData()

    for i in range(0, len(bin_heights)):
        if(centerVal[i]< 120.0):
            total = total + centerVal[i]

    for i in range(0, len(bin_heights)):
        if(centerVal[i]< 120.0):
            length = length + bin_heights[i]
    lamb = total/len(centerVal)
    print(total/len(centerVal))
    return lamb


def twob():
    fit = True
    lamb = twoa()
    print(lamb)
    totalarea = 0
    estimatedArea = 0
    AVar=1
    step = 1
    vals,centerVal,bin_heights,bin_edges = loadData()
    massVal = np.arange(104.0,155.0,step)
    print(massVal)
    for i in range(0,len(bin_heights)):
        totalarea = bin_heights[i]*1.7 + totalarea
    
    print(totalarea)
    while(fit):
        print(estimatedArea,"Estimated Area")
        print(totalarea,"totalarea")
        estimatedArea = 0
        print("-----------------------")
        for i in range(0,len(massVal)):
            estimatedArea = AVar*np.exp(-massVal[i]/lamb)*step+estimatedArea

        if(estimatedArea<totalarea):
            AVar = AVar +1
        else:
            AVar = AVar - 1
            fit = False

    std = np.sqrt(bin_heights)

    g = plt.figure(2)
    plt.plot(massVal,AVar*np.exp(-massVal/lamb))
    plt.hist(vals, range = [104.0,155.0], bins = 30)
    plt.errorbar(centerVal, bin_heights, yerr=std, fmt='r.')
    g.show()
    input()

def twob_():
    fit = True
    lamb = twoa()
    print(lamb)
    totalarea = 0
    estimatedArea = 0
    AVar=1
    step = 1
    vals,centerVal,bin_heights,bin_edges = loadData()
    massVal = np.arange(104.0,155.0,step)
    print(massVal)
    for i in range(0,len(bin_heights)):
        totalarea = bin_heights[i]*1.7 + totalarea
    
    print(totalarea)
    while(fit):
        print(estimatedArea,"Estimated Area")
        print(totalarea,"totalarea")
        estimatedArea = 0
        print("-----------------------")
        for i in range(0,len(massVal)):
            estimatedArea = AVar*np.exp(-massVal[i]/lamb)*step+estimatedArea

        if(estimatedArea<totalarea):
            AVar = AVar +1
        else:
            AVar = AVar - 1
            fit = False
    
    print(AVar, "A")

    return lamb,AVar,vals,centerVal,bin_heights,bin_edges

def twoc():
    lamb,AVar,vals,centerVal,bin_heights,bin_edges = twob_()
    rangeVal = [104.0,155.0]
    para = sht.get_B_chi(vals,rangeVal,30,AVar,lamb)
    print(para)

def twod():
    lamb,AVar,vals,centerVal,bin_heights,bin_edges = twob_()
    std = np.sqrt(bin_heights)
    step = 1
    massVal = np.arange(104.0,155.0,step)
    chiArray = [9999,0,0]
    chi = 0
    A = AVar 
    if(A < 0):
        A = 1
    lamb = lamb - 10
    for i in range(0,50000):
        chi = 0
        Av = A + i
        if(i%10000 == 0):
            print(i)
        for i in range(0,20):
            lambv = lamb +i
            expected = Av*np.exp(-bin_edges/lambv)
            for i in range(0,len(bin_heights)):
                chi = (bin_heights[i]- expected[i])**2/expected[i] + chi
            chi = chi/float(30-2)
            if(chiArray[0]> chi):
                chiArray = [chi,Av,lambv]
    print(chiArray)
    o = plt.figure(2)
    plt.plot(massVal,chiArray[1]*np.exp(-massVal/chiArray[2]))
    #plt.plot(massVal,61952*np.exp(-massVal/29.480))
    plt.hist(vals, range = [104.0,155.0], bins = 30)
    plt.errorbar(centerVal, bin_heights, yerr=std, fmt='r.')
    o.show()
    input()
   
def fit_curve(x,A,lamb):
    curve = A*np.exp(-x/lamb)
    return curve

def curve_fit():
    vals,centerVal,bin_heights,bin_edges = loadData()
    std = np.sqrt(bin_heights)
    print(len(centerVal))
    print(len(bin_heights))
    initial_guess = [28940,35]
    po1,cov1 = sp.optimize.curve_fit(fit_curve,centerVal,bin_heights,initial_guess,std)
    print(po1)

twod()
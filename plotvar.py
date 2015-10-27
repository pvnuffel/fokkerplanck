# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:50:43 2015

@author: pieter
"""

import numpy as np
import matplotlib.pylab as plt
import sys, glob

np.savetxt('results/data/Jv.out', Jv_norm)
Jv_norm = np.loadtxt('results/data/Jv.out')

M=  len( Jv_norm)
Jv_var = scipy.zeros(M)
average= Jv_norm[-1]
for m in range(1,M):
    a= np.power((Jv_norm[m] - average),2) 
    Jv_var[m]= np.sqrt(a)

ax=plt.subplot(111)
ax.set_xlim(1, 5000)
ax.set_ylim(0, np.max(Jv_var))

plot_var = plt.plot(Jv_var)
plt.savefig("results/plots/variance_Jv.pdf")

plt.show()

"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np
import matplotlib.pylab as plt

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


if __name__=="__main__":

    residual_directsim = np.loadtxt('Newton/residual_sim_N1e6.out')/sqrt(340) 
    residual1 = np.loadtxt('Newton/residual1_N1e6.out') 
    residual2 = np.loadtxt('Newton/residual2_N1e6.out') 
    residual3 = np.loadtxt('Newton/residual3_N1e6.out')  
    
#    
#
#            
#        residual_directsim[Dti] = norm( rho_ss- rho_Dt)
#        residual2[Dti] = norm( States[2] - rho_ss)/sqrt(len(States[2] - rho_ss))
#        residual1[Dti] = norm(States[1]-rho_ss)/sqrt(len(States[1] - rho_ss))
#        residual3[Dti] = norm(States[3]-rho_ss)/sqrt(len(States[3] - rho_ss))
  
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(Dtlist, residual1, "g--")
    plt.plot(Dtlist, residual2, "r--" )
    plt.plot(Dtlist, residual3, "y--" )
    plt.plot(Dtlist, residual_directsim, "c--" )
    plt.xlabel(r'$\Delta T$', fontsize = 14)
    plt.ylabel(r'$||  \rho^{*}   - \frac{ \exp{\left[-2\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}  ||$'  , fontsize = 15)
    plt.savefig("plots/Newton_sde_res(Dt)_N1e6_dt=1e-4.pdf")
              

        

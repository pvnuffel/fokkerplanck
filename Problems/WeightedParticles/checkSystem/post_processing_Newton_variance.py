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
from matplotlib.path import Path #to draw trriangles indicating order of convergence
import matplotlib.patches as patches #to draw trriangles indicating order of convergence

from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt


if __name__=="__main__":
   

    
    Nlist = scipy.array([1e3, 1e4, 1e5, 1e6])
    Nlist_inv = 1.0/Nlist
    rho_sq_tol2 = np.loadtxt('Newton/rho_sq.out')         #M=50
    rho_sq_tol4 = np.loadtxt('Newton/rho_sq_tol_e-4.out')  
    rho_sq_tol3 = np.loadtxt('Newton/rho_sq_tol_e-3.out') 
    rho_sq_tol5 = np.loadtxt('Newton/rho_sq_tol_e-5.out')  
    #np.loadtxt('25-11-Rho_pde_dx_e-4_eps-6.out') ,  np.loadtxt('25-11-Rho_pde_dx_e-4_eps-7.out') ]  #these data should be the same
#    Jv_pde= np.loadtxt('data/25-11-Jv_pde_dx_e-4_eps-5.out') 
   # Jv_pde_dxbig= np.loadtxt('data/Jv_pde_dx_5e-4_vsin.out') 
    E_rho_tol2 = np.loadtxt('Newton/E_rho.out')       #M=50
    E_rho_tol4 = np.loadtxt('Newton/E_rho_tol_e-4.out')  
    E_rho_tol3 = np.loadtxt('Newton/E_rho_tol_e-3.out') 
    E_rho_tol5 = np.loadtxt('Newton/E_rho_tol_e-5.out') 
  #  Jv_pde = Jv_pde_list[1]
    
    #SDE       
    variance =  scipy.zeros(len(Nlist))
    variance4 =  scipy.zeros(len(Nlist))
    variance3 =  scipy.zeros(len(Nlist))
    variance1 =  scipy.zeros(len(Nlist))
    variance5 =scipy.zeros(len(Nlist))
    
    for n in range (0, len(Nlist)):
        variance[n] = rho_sq_tol2[n] - norm(E_rho_tol2[n])**2
        variance4[n] = rho_sq_tol4[n] - norm(E_rho_tol4[n])**2
        variance3[n] =  rho_sq_tol3[n] - norm(E_rho_tol3[n])**2 
        variance1[n] =  rho_sq_tol1[n] - norm(E_rho_tol1[n])**2 
        variance5[n] =  rho_sq_tol5[n] - norm(E_rho_tol5[n])**2 
    
    log_flag = 'log'
  #  log_flag = 'linear'
    save_flag = True

    bins = len(E_rho[-1])
    
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])   
    plot_var = plt.plot(Nlist_inv, variance/bins, label=r'$\epsilon_{GMRES}=10^{-2}$', linewidth=2  )
    plot_var = plt.plot(Nlist_inv, variance3/bins,  label=r'$\epsilon_{GMRES}=10^{-3}$',linewidth=2  )
    plot_var = plt.plot(Nlist_inv, variance4/bins,  label=r'$\epsilon_{GMRES}=10^{-4}$',linewidth=2  )
    plot_var = plt.plot(Nlist_inv, variance5/bins,  label=r'$\epsilon_{GMRES}=10^{-5}$',linewidth=2  )
   # plot_var = plt.plot(Nlist_inv, variance1/bins, linewidth=2  )
    plt.ylabel(r'Var ($\rho^*_{k=3}$)', fontsize = 16)
    plt.xlabel('$1/N$', fontsize = 16) 
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.legend([plot_var], loc='best')
    plt.legend(bbox_to_anchor=(0.4,1),  numpoints = 1 )
    if(save_flag): plt.savefig("plots/Variance_3th_newton_it_Dt1e-1.pdf")
    plt.show()

            

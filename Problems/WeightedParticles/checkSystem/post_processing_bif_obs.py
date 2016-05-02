# -*- coding: utf-8 -*-
import scipy
import numpy as np
from pylab import *
from scipy.linalg import norm
"""
Created on Wed Mar 23 20:01:01 2016

@author: pieter
"""

dx=1e-2
xL = -1.7
xR = 1.7
grid = scipy.arange(xL+dx/2.,xR,dx)

M=10
rho_mean = np.loadtxt('bifurcation_data/fixed_states_sde_Ne6_mean_M10.txt')
y_mean =  np.loadtxt('bifurcation_data/norm_fixed_states_sde_Ne6_mean_M10_RL.txt')
y_var = np.loadtxt('bifurcation_data/variance_fixed_states_sde_Ne6_mean_M10_RL.txt')
y_error = sqrt(y_var)/sqrt(M)

#STOCHASTIC

save_flag= False
N_points = len(y_mean)

sigma_list= scipy.arange(2, 0.9, -0.1)
        
for i in range (N_points):
    par_label = r'$\sigma= %g $' % sigma_list[i]
    plt.plot(grid, rho_mean[i], label=par_label)   
plt.xlabel(r'$\sigma$', fontsize = 15)
plt.legend()
if save_flag:
    plt.savefig('plots/bif/fixed_states_sde(sigma)_Ne6_mean_M10.pdf')  
plt.show()

#ANALYTIC
mu=1.0
sigma_list= scipy.arange(1.5, 0.4, -0.1)
sigma_list= scipy.arange(2, 0.9, -0.1)
D_list= scipy.arange(1, 0, -0.1)
y_anal = []
rho_norm_D = []
for i in range ( len(sigma_list)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2 )*2*mu/(sigma_list[i]**2)) /norm_c
    y_anal.append(norm(rho_ss))

sigma_list_fine= scipy.arange(2.1, 0.9, -0.001)
y_anal_fine = []
for i in range ( len(sigma_list_fine)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list_fine[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2 )*2*mu/(sigma_list_fine[i]**2)) /norm_c
    y_anal_fine.append(norm(rho_ss))



#PLOT BIFURCATION 

    
plt.plot( sigma_list_fine,y_anal_fine,  'r', markersize=4,label= 'Analytic solution')
plt.errorbar(sigma_list,  y_mean, yerr=y_error,   fmt='bo', ecolor='b', markersize=3, label='Newton-Krylov solution')
plt.ylabel(r'$\left\||\rho^*\right\||$', fontsize = 15)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.xlim([0.9,2.1])
plt.legend(prop={'size':8})

#plt.plot(x_sde, y_sde, 'go', markersize=4)

if save_flag:
    plt.savefig('plots/bif/bifurcation_sde_Ne6_anal(sigma)_RL.pdf')

plt.show()


plt.xlim([0.9,2.1])
plt.ylabel(r'$\left\||\rho^*\right\||$', fontsize = 10)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.errorbar(sigma_list, y_mean-y_anal, yerr=y_error,   fmt='bo', ecolor='b', markersize=3) # Mean of M realizations

if save_flag:
    plt.savefig('plots/bif/bifurcation_bias(sigma)_sde_Ne6_mean_M10_RL.pdf')

plt.show()
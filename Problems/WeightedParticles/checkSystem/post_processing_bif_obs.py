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
y_mean =  np.loadtxt('bifurcation_data/norm_fixed_states_sde_Ne6_mean_M10_LR.txt')
y_var = np.loadtxt('bifurcation_data/variance_fixed_states_sde_Ne6_mean_M10_LR.txt')
y_error = sqrt(y_var)/sqrt(M)

#STOCHASTIC

save_flag= True
N_points = len(y_mean)

#sigma_list= scipy.arange(2, 0.9, -0.1)   #voor RL
sigma_list= scipy.arange(1.0, 2.1, 0.1)  #voor LR
        
for i in range (N_points):
    if (i % 2 == 0): 
        par_label = r'$\sigma= %g $' % sigma_list[i]
        plt.plot(grid, rho_mean[i], label=par_label)   
plt.xlim([-1.7,1.7])
plt.ylim([0,0.6])
plt.xlabel('$x$', fontsize = 12)
plt.ylabel(r'$\rho^*$', fontsize = 12)
plt.legend(bbox_to_anchor=(1,0.6),  numpoints = 1 )
plt.legend(prop={'size':9})
if save_flag:
    plt.savefig('plots/bif/fixed_states_sde(sigma)_Ne6_mean_M10_LR.pdf')  
plt.show()

#ANALYTIC
mu=1.0
#sigma_list= scipy.arange(1.5, 0.4, -0.1)
#sigma_list= scipy.arange(2, 0.9, -0.1)
D_list= scipy.arange(1, 0, -0.1)
y_anal = []
rho_norm_D = []
for i in range ( len(sigma_list)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2 )*2*mu/(sigma_list[i]**2)) /norm_c
    y_anal.append(norm(rho_ss))
    if (i % 2 == 0): 
        par_label = r'$\sigma= %g $' % sigma_list[i]
        plt.plot(grid, rho_ss, label=par_label)

plt.xlim([-1.7,1.7])
plt.ylim([0,0.6])
plt.xlabel('$x$', fontsize = 12)
plt.ylabel(r'$\rho^*$', fontsize = 12)
plt.legend(bbox_to_anchor=(1,0.6),  numpoints = 1 )
plt.legend(prop={'size':9})
if save_flag:
    plt.savefig('plots/bif/fixed_states(sigma)_analytic.pdf')  
plt.show()


sigma_list_fine= scipy.arange(2.1, 0.9, -0.001)
y_anal_fine = []
for i in range ( len(sigma_list_fine)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list_fine[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2 )*2*mu/(sigma_list_fine[i]**2)) /norm_c
    y_anal_fine.append(norm(rho_ss))




#PLOT BIFURCATION 

ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])       
plt.plot( sigma_list_fine,y_anal_fine,  'r', markersize=4,label= 'Analytic solution')
plt.errorbar(sigma_list,  y_mean, yerr=y_error,   fmt='bo', ecolor='b', markersize=3, label='Newton-Krylov solution')
plt.ylabel(r'$\left\||\rho^*\right\||_2$', fontsize = 12)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.xlim([0.9,2.1])
plt.legend(prop={'size':8})

#plt.plot(x_sde, y_sde, 'go', markersize=4)

if save_flag:
    plt.savefig('plots/bif/bifurcation_sde_Ne6_anal(sigma)_LR.pdf')

plt.show()


plt.xlim([0.9,2.1])
plt.ylabel(r'Bias($\left\||\rho^*\right\||_2)$', fontsize = 12)
plt.xlabel(r'$\sigma$', fontsize = 12)
#plt.yscale('log')
plt.errorbar(sigma_list, y_mean-y_anal, yerr=y_error,   fmt='bo', ecolor='b', markersize=3) # Mean of M realizations

if save_flag:
    plt.savefig('plots/bif/bifurcation_bias(sigma)_sde_Ne6_mean_M10_LR.pdf')

plt.show()
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

#STOCHASTIC

save_flag= False
branch = branch_list[-1]

bif_par = scipy.zeros( len(branch.points))

for i in range ( len(branch.points)):
    bif_par[i] = branch.points[i].lambd[1]
    
for i in range (len(branch.points)):
  #  par_label = r'$a = %g $' % bif_par[i]
    par_label = r'$\sigma= %g $' % bif_par[i]
    plt.plot(grid, branch.points[i].u,  label=par_label)
plt.xlabel(r'$\sigma$', fontsize = 15)
plt.legend()
if save_flag:
    plt.savefig('plots/bif/fixed_states_sde(sigma)_Ne5_2.pdf')
    
plt.show()


#for b in range ( len(branch_list)):
#    plt.plot(grid, branch_list[b].points[-1].u, label=b)
#plt.xlabel(r'$\sigma$', fontsize = 15)
#plt.show()
#    
    
#x, y=  branch.getPlot()

# MEAN FIXED STATE
M= len(branch_list)
rho_mean =  scipy.zeros((len(branch.points), len(branch.points[-1].u)))

for b in range (M):
    for i in range (len(branch_list[b].points)):
        rho_mean[i] =  rho_mean[i] + branch_list[b].points[i].u/M

if save_flag:
    np.savetxt('bifurcation_data/fixed_states_sde_Ne5_mean_M10_2.txt',rho_mean)
        
        
for i in range (len(rho_mean)):
    par_label = r'$\sigma= %g $' % bif_par[i]
    plt.plot(grid, rho_mean[i], label=par_label)   
plt.xlabel(r'$\sigma$', fontsize = 15)
plt.legend()
if save_flag:
    plt.savefig('plots/bif/fixed_states_sde(sigma)_Ne5_mean_M10_2.pdf')  
plt.show()

#ANALYTIC

def resize( original_vector, new_size):
    resize_factor = int (len(original_vector)/new_size) 
    #print "Discretisation for solving sde is ",  resize_factor , " times coarser than the discretisation for solving the pde"  
    new_vector = scipy.zeros(new_size)
    for i in range (0,new_size):
        bin_av = 0
        for j in range (0,  resize_factor ):
            bin_av = bin_av+ original_vector[i*resize_factor+j]
            new_vector[i] = bin_av/resize_factor
    return new_vector
    






sigma_list_fine= scipy.arange(1, 2.1, 0.001)
y_anal_fine = []
for i in range ( len(sigma_list_fine)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list_fine[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2 )*2*mu/(sigma_list_fine[i]**2)) /norm_c
    y_anal_fine.append(norm(rho_ss))



mu=1.0
sigma_list= scipy.arange(1.5, 0.4, -0.1)
#sigma_list= scipy.arange(2, 0.9, -0.1)
sigma_list= scipy.arange(1, 2.1, 0.1)
D_list= scipy.arange(1, 0, -0.1)

rho_norm_D = []
y_anal_test =[] #without rescaling

for i in range ( len(sigma_list)):
    norm_c = sum( np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list[i]**2)))*dx
    rho_ss = np.exp( (-grid**4 + grid**2)*2*mu/(sigma_list[i]**2))/norm_c
    y_anal_test.append(norm(rho_ss))     

dx=1e-6
y_anal = []
grid_fine = scipy.arange(xL+dx/2.,xR,dx)
for i in range ( len(sigma_list)):
    norm_c = sum( np.exp( (-grid_fine**4 + grid_fine**2)*2*mu/(sigma_list[i]**2)))*dx
    rho_ss_fine = np.exp( (-grid_fine**4 + grid_fine**2)*2*mu/(sigma_list[i]**2))/norm_c
    rho_ss_coarse = resize(rho_ss_fine, len(grid))
    y_anal.append(norm(rho_ss_coarse))               #nauwkeuriger dan rho_ss
     

    #plt.plot(rho_ss)

#for i in range ( len(D_list)):
#    norm_c = sum( np.exp( (-grid**4 + grid**2)*mu/D_list[i]))*dx
#    rho_ss = np.exp( (-grid**4 + grid**2 )*mu/D_list[i]) /norm_c
#    rho_norm_D.append(norm(rho_ss))
#    plt.plot(rho_ss)
#plt.show()

#for i in range (len(rho_norm_D)):
#    plt.plot( D_list,rho_norm_D, 'o')
#    plt.xlabel(r'$D$', fontsize = 12)
    
  
    
##DETERMINISTIC
#bif_par_pde = zeros( len(branch_pde.points))
#
#for i in range ( len(branch_pde.points)):
#    bif_par_pde[i] = branch_pde.points[i].lambd[1] 
#    
#for i in range (len(branch_pde.points)):
#  #  par_label = r'$\alpha = %g $' % bif_par_pde[i]
#    par_label = r'$D = %g $' % bif_par_pde[i]
#    plt.plot(branch_pde.points[i].u, label=par_label)
#
#plt.legend()
#if save_flag:
#    plt.savefig('plots/fixed_states_pde(D).pdf')
#plt.show()
#




#PLOT BIFURCATION 


new_branch_list=branch_list
for b in range (M):
    print 'branch ', b, ' has ' , len(branch_list[b].points), ' points'
    if (len(new_branch_list[b].points) != len(branch.points)):
        new_branch_list.remove(new_branch_list[b])
        print  'branch ', b, ' has been removed'
        
y_s1 =  scipy.zeros(len(branch.points))
y_s2 =  scipy.zeros(len(branch.points))
x_sde= branch_list[-1].getPlot()[0]

y_mean =  scipy.zeros(len(branch.points))

for b in range (M):
    y_sde=branch_list[b].getPlot()[1]
    y_s1 += y_sde
    y_s2 += y_sde*y_sde 
    plt.plot(x_sde, y_sde, 'go', markersize=4)   #plot all branches (all realizations)
    
plt.show()
 #   par_label = r'$D= %g $' % bif_par[b]
#    plt.plot(branch.points[b].u, label=par_label)
    

y_mean = y_s1/M
y_var = y_s2/M - y_mean*y_mean
y_sigma= sqrt(y_var)
y_error = y_sigma/sqrt(M)

if save_flag:
    np.savetxt('bifurcation_data/norm_fixed_states_sde_Ne5_mean_M10_LR2.txt',y_mean)
    np.savetxt('bifurcation_data/variance_fixed_states_sde_Ne5mean_M10_LR2.txt',y_var)


   
    
plt.plot( sigma_list_fine,y_anal_fine,  'r', markersize=4,label= 'Analytic solution')
plt.errorbar(x_sde, y_mean, yerr=y_error,   fmt='bo', ecolor='b', markersize=3, label='Newton-Krylov solution')
plt.ylabel(r'$\left\||\rho^*\right\||$', fontsize = 15)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.xlim([0.9,2.1])

#plt.plot(x_sde, y_sde, 'go', markersize=4)

if save_flag:
    plt.savefig('plots/bif/bifurcation_sde_Ne5_anal(sigma)_LR2.pdf')

plt.show()

plt.errorbar(x_sde, y_sde-y_anal, yerr=y_sigma,   fmt='bo', ecolor='b', markersize=3) # For 1 realization
plt.xlim([0.9,2.1])

plt.ylabel(r'Bias($\left\||\rho^*\right\||_2$)', fontsize = 12)
plt.xlabel(r'$\sigma$', fontsize = 11)

if save_flag:
    plt.savefig('plots/bif/bifurcation_bias(sigma)_sde_Ne5_LR2.pdf')
plt.show()


plt.xlim([0.9,2.1])
plt.ylabel(r'Bias($\left\||\rho^*\right\||_2)$', fontsize = 12)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.errorbar(x_sde, y_mean-y_anal, yerr=y_error,   fmt='bo', ecolor='b', markersize=3) # Mean of M realizations

if save_flag:
    plt.savefig('plots/bif/bifurcation_bias(sigma)_sde_Ne5_mean_M10-1_LR2.pdf')

plt.show()

plt.xlim([0.9,2.1])
plt.ylabel('Discretization error', fontsize = 12)
plt.xlabel(r'$\sigma$', fontsize = 11)
plt.plot(sigma_list, np.array(y_anal_test)-np.array(y_anal), 'o', markersize=3) # Mean of M realizations

plt.show()
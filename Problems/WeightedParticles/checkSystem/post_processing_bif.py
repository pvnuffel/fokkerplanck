# -*- coding: utf-8 -*-
import scipy
"""
Created on Wed Mar 23 20:01:01 2016

@author: pieter
"""

#prev_point = branch.points[0].u
#last_point = branch.points[1].u

#STOCHASTIC

save_flag= False
bif_par = zeros( len(branch.points))

for i in range ( len(branch.points)):
    bif_par[i] = branch.points[i].lambd[0]
    
for i in range (len(branch.points)):
    par_label = r'$a = %g $' % bif_par[i]
  #  par_label = r'$D = %g $' % bif_par[i]
    plt.plot(branch.points[i].u, label=par_label)
    
plt.legend()
if save_flag:
    plt.savefig('plots/fixed_states_sde(a)_1.pdf')
    
plt.show()
    
#x, y=  branch.getPlot()


plt.show()
    
#DETERMINISTIC
bif_par_pde = zeros( len(branch_pde.points))

for i in range ( len(branch_pde.points)):
    bif_par_pde[i] = branch_pde.points[i].lambd[0]
    
for i in range (len(branch_pde.points)):
  #  par_label = r'$\alpha = %g $' % bif_par_pde[i]
    par_label = r'$a = %g $' % bif_par_pde[i]
    plt.plot(branch_pde.points[i].u, label=par_label)

plt.legend()
if save_flag:
    plt.savefig('plots/fixed_states_pde(a)_1.pdf')
plt.show()
    
#x_pde, y_pde=  branch_pde.getPlot()
plt.plot(x_pde, y_pde, 'ro') #, markersize=3)
plt.ylabel(r'$\left\||\rho^*\right\||$', fontsize = 16)
plt.xlabel(r'$a$', fontsize = 12)
plt.savefig('plots/bifurcation_pde(D).pdf')

plt.plot(x_sde, y_sde, 'bo')

if save_flag:
    plt.savefig('plots/bifurcation_pde_1.pdf')

plt.show()

plt.plot(x_sde, abs(y_sde-y_pde), 'bo')

if save_flag:
    plt.savefig('plots/bifurcation_bias_1.pdf')

plt.show()
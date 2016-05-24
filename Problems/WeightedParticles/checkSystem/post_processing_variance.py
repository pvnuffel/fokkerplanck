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
plt.rc('text', usetex=True)

if __name__=="__main__":
    D = 1./2
    Dt = 1e-2    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
   # param['eps'] = 1e-5
    #param['Dt'] = Dt
    #param['dt']=dt         

    
    #SDE       
#    rho_sq = np.loadtxt('data/obs/25-11-rho_sq_eps-5.out')
#    E_rho = np.loadtxt('data/obs/25-11-E_rho_eps-5.out')
#    sq_E_rho  =np.loadtxt('data/obs/25-11-sq_E_rho_eps-5.out')
#    
    rho_sq = np.loadtxt('data/09-05-rho_sq_N=e^3-7.out')
    E_rho = np.loadtxt('data/09-05-E_rho_N=e^3-7.out')
    #sq_E_rho  =np.loadtxt('data/09-05-sq_E_rho.out')
    sq_E_rho  =np.loadtxt('data/09-05-sq_E_rho_N=e^3-7.out')
  
    
#            np.savetxt('data/09-05-rho_sq.out' , rho_sq)
#        np.savetxt('data/09-05-E_rho.out' , E_rho )
#        np.savetxt('data/09-05-sq_E_rho.out',  sq_E_rho )
        
     
    Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000])
    Nlist = scipy.array([1e3, 1e4, 1e5, 1e6, 1e7])
    Nlist_inv = 1.0/Nlist
    
    log_flag = 'log'
  #  log_flag = 'linear'
    save_flag = False
 

    bins = len(E_rho[-1])
    
    variance = (rho_sq - sq_E_rho)
    sdev = sqrt(variance)
    
    points_o1var = [[1e-6, 2e-2], [1e-6, 1e-1], [5e-6, 1e-1]] 
 #   points_o1var = [[1e-5, 1e-2], [2e-5, 2e-2], [1e-5, 2e-2]] 
    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
        
    plot_var = plt.plot(Nlist_inv, variance, label = '$Variance$', linewidth=2)
    plt.ylabel(r'Var$(\Phi^N_T(\rho))$' , fontsize = 12)
    plt.xlabel('$1/N$', fontsize = 12)
    plt.xscale(log_flag)
    plt.yscale(log_flag)
    plt.gca().add_patch(triangle_o1var)
    order= r'$\mathcal{O}(1/N)$'
    plt.annotate(order,  xy=(2e-6, 2.5e-2), xytext=(1.1e-6, 5.7e-2), fontsize=11, color='grey')
    plt.savefig('plots/Variance_on_cts(N-1).pdf')
    plt.show()
    

#    plot_mse = plt.plot(Nlist, (sde_rho_sq -pde_rho_sq)/bins, label = '$RMS^2$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'MSE( $\rho$)' , fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
# 
#    plot_bias =  plt.plot(Nlist, (abs(pde_rho_sq - sq_E_rho))/bins, label = '$Bias$')
#    plt.xlabel('$N$', fontsize = 16)
#    plt.ylabel(r'Bias $(\rho)$', fontsize = 16)
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.show()
    
 #   plot_rho = plt.plot( rho_coarse, label = '$\$')
 #   plot_rho =  plt.plot(rho_Dt_sde, label = '$\$')
  #  plt.ylabel(r' $\rho$', fontsize = 16)
 #   plt.xlabel('$x$', fontsize = 16)
  #  plt.show()

#    lines = ["-","--","-.",":"]    
#    
#    points_o1var = [[1e-4, 2e-2], [1e-4, 1e-1], [5e-4, 1e-1]] 
# #   points_o1var = [[1e-5, 1e-2], [2e-5, 2e-2], [1e-5, 2e-2]] 
#    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
#    
#    
#    #VariancePlot
#       
#    for eps_i in range (0,len(eps_list)):
#        plot_var = plt.plot(Nlist_inv, (sde_Jv_sq[eps_i] -sq_E_Jv[eps_i])/bins, lines[eps_i],
#                            label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] , linewidth=3)
#    plt.ylabel('Var ($\mathbf{\hat{Jv}} $)', fontsize = 16)
#    plt.xlabel('$1/N$', fontsize = 16) 
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    plt.legend([plot_var], loc='best')
#    plt.legend(bbox_to_anchor=(1, 0.5), numpoints = 1 )
#    plt.gca().add_patch(triangle_o1var)
#    order= r'$\mathcal{O}(1/N)$'
#    plt.annotate(order,  xy=(2e-4, 2e-2), xytext=(1.1e-4, 5.2e-2), fontsize=11, color='grey')
#    if(save_flag): plt.savefig("plots/Var_N_eps.pdf")
#    plt.show()
#    

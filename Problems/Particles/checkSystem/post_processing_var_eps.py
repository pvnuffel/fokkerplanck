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
    eps_list_exponents=[4,5, 6,7]
    eps_list= [1e-4, 1e-5, 1e-6, 1e-7]


    #SDE       
   # sde_rho_sq = np.loadtxt('data/25-11-rho_sq_eps-5.out')
   # E_rho = np.loadtxt('data/25-11-E_rho_eps-5.out')
   # sq_E_rho  =np.loadtxt('data/25-11-sq_E_rho_eps-5.out')
    
    sde_Jv_sq = [ np.loadtxt('nw-Jv_sq_eps-4.out'), np.loadtxt('nw-Jv_sq_eps-5.out'),
           np.loadtxt('nw-Jv_sq_eps-6.out'), np.loadtxt('nw-Jv_sq_eps-7.out') ] 
#    E_Jv = [np.loadtxt('data/eps_fin-E_Jv_eps-4.out'),  np.loadtxt('data/eps_fin-E_Jv_eps-5.out'), 
#            np.loadtxt('data/eps_fin-E_Jv_eps-6.out'), np.loadtxt('data/eps_fin-E_Jv_eps-7.out') ]    # [eps_i x N_i x dx] = [4 x 7 x 3400]
    sq_E_Jv =  [np.loadtxt('nw-sq_E_Jv_eps-4.out') , np.loadtxt('nw-sq_E_Jv_eps-5.out' ) ,
                 np.loadtxt('nw-sq_E_Jv_eps-6.out' ),np.loadtxt('nw-sq_E_Jv_eps-7.out' )]
                 
                 
    sde_Jv_sq_weighted = [ np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-Jv_sq_eps-4.out'), np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-Jv_sq_eps-5.out'), 
                         np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-Jv_sq_eps-6.out'), np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-Jv_sq_eps-7.out') ] 
    sq_E_Jv_weighted =  [np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-sq_E_Jv_eps-4.out' ), np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-sq_E_Jv_eps-5.out' ), 
                np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-sq_E_Jv_eps-6.out' ),np.loadtxt('../../WeightedParticles/checkSystem/data/eps_fin-sq_E_Jv_eps-7.out' )]
                

 
    #rho_Dt_sde = E_rho[-1]
   # Jv_sde = E_Jv[-1][-1]  #dim n_x

 
    Nlist = scipy.array([1000, 2000,4000,8000,16000,32000,64000, 128000, 256000])
    Nlist_inv = 1.0/Nlist
    
    log_flag = 'log'
  #  log_flag = 'linear'v\
    save_flag = True
    bins = 340

    lines = ["g-","b--","r-.","c:"]    

   # lines = ["-","-","-","-"]   
    
    points_o1var = [[1e-5, 1e5], [1e-4, 1e5], [1e-4, 1e6]] 
 #   points_o1var = [[1e-5, 1e-2], [2e-5, 2e-2], [1e-5, 2e-2]] 
    triangle_o1var = plt.Polygon(points_o1var, fill=None  ,edgecolor='grey') 
    
    
    #VariancePlot
       
    for eps_i in range (0,len(eps_list)):
        plot_var = plt.plot(Nlist_inv, (sde_Jv_sq[eps_i] -sq_E_Jv[eps_i])/bins,  lines[eps_i],
                            label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] , linewidth=3)
        plot_weighted = plt.plot(Nlist_inv, (sde_Jv_sq_weighted[eps_i] -sq_E_Jv_weighted[eps_i])/bins, lines[eps_i], linewidth=3.5)
    plt.ylabel('Var ($\mathbf{\hat{Jv}} $)', fontsize = 16)
    plt.xlabel('$1/N$', fontsize = 16) 
    plt.xscale(log_flag)
    plt.yscale(log_flag)
   # plt.legend([plot_var], loc='best')
    #plt.legend(bbox_to_anchor=(0.95, 1), numpoints = 1 )
    plt.gca().set_ylim([1e-4, 1e15])   
    plt.gca().add_patch(triangle_o1var)
    plt.legend([plot_var], loc='best')
    plt.legend(bbox_to_anchor=(1, 0.65), numpoints = 1 )    
    order= r'$\mathcal{O}(1/N)$'
    plt.annotate(order,  xy=(4.8e-5, 1.5e5), xytext=(4.8e-5, 1.5e5), fontsize=9, color='grey')
    plt.annotate('with variance reduction',  xy=(4e-6, 1e-1), xytext=(4e-6, 1e-1), fontsize=12,rotation=5.5, color='grey')
    if(save_flag): plt.savefig("plots/Var_N_eps_nw.pdf")
    plt.show()
    
#    for eps_i in range (0,len(eps_list)):
#        plot_var = plt.plot(Nlist_inv, (rho_sq -sq_E_rho)/bins, lines[eps_i],
#                            label =r'$\varepsilon=10^{-%d}$' %eps_list_exponents[eps_i] , linewidth=2)
#    plt.ylabel('Var ($\mathbf{\hat{rho}} $)', fontsize = 16)
#    plt.xlabel('$1/N$', fontsize = 16) 4.4e-5, 1.5e5
#    plt.xscale(log_flag)
#    plt.yscale(log_flag)
#    
#    
    

    
   # plt.legend([plot_var], loc='best')
    #plt.legend(bbox_to_anchor=(0.95, 1), numpoints = 1 )
  #  plt.gca().set_ylim([1e6, 1e15])   
#    plt.gca().add_patch(triangle_o1var)
  #  if(save_flag): plt.savefig("plots/Var_N_eps_nw.pdf")
    plt.show()
    
    

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

#    residual_directsim = np.loadtxt('Newton/residual_sim_N1e6.out')
#    residual1 = np.loadtxt('Newton/residual1_N1e6.out') 
#    residual2 = np.loadtxt('Newton/residual2_N1e6.out') 
#    residual3 = np.loadtxt('Newton/residual3_N1e6.out')  
#    
#    residual_directsim = np.loadtxt('Newton/residual_sim(N)_Dte-1_tol1e-3.out')
#    residual1 = np.loadtxt('Newton/residual1(N)__Dte-1_tol1e-3.out') 
#    residual2 = np.loadtxt('Newton/residual2(N)__Dte-1_tol1e-3.out') 
#    residual3 = np.loadtxt('Newton/residual3(N)_Dte-1_tol1e-3.out')  
#    

#    resnorm_t100_N1000 = np.loadtxt('Newton/resnorm_Dte-1_N1000')
#    resnorm_t100_N10000 = np.loadtxt('Newton/resnorm_Dte-1_N10000') 
#    resnorm_t100_N100000 = np.loadtxt('Newton/resnorm_Dte-1_N100000') 
#    resnorm_t100_N1000000 = np.loadtxt('Newton/resnorm_Dte-1_N1000000')  
#    resnorm_t100_N10000000= np.loadtxt('Newton/resnorm_Dte-1_N10000000')  
    
    resnorm_t10_tol2_N1000 = np.loadtxt('Newton/resnorm_Dte-2_tole-2_N1000')
    resnorm_t10_tol2_N10000 = np.loadtxt('Newton/resnorm_Dte-2_tole-2_N10000') 
    resnorm_t10_tol2_N100000 = np.loadtxt('Newton/resnorm_Dte-2_tole-2_N100000') 
    resnorm_t10_tol2_N1000000 = np.loadtxt('Newton/resnorm_Dte-2_tole-2_N1000000')  
    resnorm_t10_tol2_N10000000= np.loadtxt('Newton/resnorm_Dte-2_tole-2_N10000000')  

    
    resnorm_t10_N1000 = np.loadtxt('Newton/resnorm_Dte-2_N1000')
    resnorm_t10_N10000 = np.loadtxt('Newton/resnorm_Dte-2_N10000') 
    resnorm_t10_N100000 = np.loadtxt('Newton/resnorm_Dte-2_N100000') 
    resnorm_t10_N1000000 = np.loadtxt('Newton/resnorm_Dte-2_N1000000')  
    resnorm_t10_N10000000= np.loadtxt('Newton/resnorm_Dte-2_N10000000')  
    
    
    
#    
   # Nlist = scipy.array([1e4,5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7]) 
       
    Nlist = scipy.array([1e3, 1e4, 1e5, 1e6, 1e7])

   # fig = plt.figure()
   # ax = fig.add_subplot(1,1,1)
    ax = plt.subplot(111)  
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.82, box.height])   
          
#        residual_directsim[Dti] = norm( rho_ss- rho_Dt)
#        residual2[Dti] = norm( States[2] - rho_ss)/sqrt(len(States[2] - rho_ss))
#        residual1[Dti] = norm(States[1]-rho_ss)/sqrt(len(States[1] - rho_ss))
#        residual3[Dti] = norm(States[3]-rho_ss)/sqrt(len(States[3] - rho_ss))
  
   # plt.xscale('log')
    plt.yscale('log')
##    plt.plot(Dtlist, residual1, "g--", label='1 Newton correction')
##    plt.plot(Dtlist, residual2, "r--",label = '2 Newton corrections ' )
##    plt.plot(Dtlist, residual3, "y--" , label= '3 Newton corrections')
##    plt.plot(Dtlist, residual_directsim, "c--" , label = 'direct simulation')
#    plt.plot(Nlist, residual1, "g--", label='1 Newton correction')
#    plt.plot(Nlist, residual2, "r--",label = '2 Newton corrections ' )
#    plt.plot(Nlist, residual3, "y--" , label= '3 Newton corrections')
#    plt.plot(Nlist, residual_directsim, "c--" , label = 'direct simulation')    
#    
#  #  plt.xlabel(r'$\Delta T$', fontsize = 14)
#    plt.xlabel(r'$N$', fontsize = 14)
#    plt.ylabel(r'$||  \rho^{*}   - \frac{ \exp{\left[-2\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}  ||$'  , fontsize = 15)
#    plt.legend(prop={'size':6})
#   # plt.savefig("plots/Newton_sde_res(Dt)_N5e5dt=1e-3_tol_1e-3.pdf")
#    plt.savefig("plots/Newton_sde_res(it)_Dt=1e-1_tol_1e-3_alpha99_dxfine.pdf")   
#    
#    
      

    plt.plot(resnorm_t10_tol2_N1000,  "g-", label=r'$N=10^3$')
    plt.plot( resnorm_t10_tol2_N10000, "r-",label = r'$N=10^4$')
    plt.plot( resnorm_t10_tol2_N100000,"y-" , label= r'$N=10^5$')
    plt.plot(  resnorm_t10_tol2_N1000000, "c-" , label = r'$N=10^6$')    
    plt.plot(  resnorm_t10_tol2_N10000000, "b-" , label = r'$N=10^7$')   
    
    
    plt.plot( resnorm_t10_N1000,  "g--")
    plt.plot( resnorm_t10_N10000, "r--")
    plt.plot( resnorm_t10_N100000,"y--" )
    plt.plot( resnorm_t10_N1000000, "c--"  )    
    plt.plot(  resnorm_t10_N10000000, "b--" )   
    
  #  plt.xlabel(r'$\Delta T$', fontsize = 14)
    plt.xlabel(r'$k$', fontsize = 14)
    plt.ylabel(r'$||  \rho - \Phi^N_T(\rho)||$'  , fontsize = 15)

  #  plt.title(r'$\Delta T = 0.01 , tol(-) = 10^{-2}, tol(--)=10^{-3} $')
    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2,3,4)
    #Create custom artists
    tol3Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
    tol2Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='-')
    
    ax.legend( numpoints = 1,  prop={'size':8} )
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[tol3Artist,tol2Artist],
          [label for i,label in enumerate(labels) if i in display]+[r'$\epsilon_{GMRES} = 10^{-3}$', r'$\epsilon_{GMRES} = 10^{-2}$'], bbox_to_anchor=(1.39
          , 1), prop={'size':10})
  
   # ax.legend(prop={'size':6})
   # plt.savefig("plots/Newton_sde_res(Dt)_N5e5dt=1e-3_tol_1e-3.pdf")
    plt.savefig("plots/Newton_sde_res(k)_Dt_e-2_alpha99.pdf")   

        


        

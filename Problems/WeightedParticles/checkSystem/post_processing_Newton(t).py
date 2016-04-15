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

    residual_directsim = np.loadtxt('Newton/residual_sim(t)_Ne5_tol1e-2.out')
    residual1 = np.loadtxt('Newton/residual1(t)_Ne5_tol1e-2.out') 
    residual2 = np.loadtxt('Newton/residual2(t)_Ne5_tol1e-2.out') 
    residual3 = np.loadtxt('Newton/residual3(t)_Ne5_tol1e-2.out') 
    
    residual_directsim_Ne6 = np.loadtxt('Newton/residual_sim(t)_Ne6_tol1e-2.out')
    residual1_Ne6  = np.loadtxt('Newton/residual1(t)_Ne6_tol1e-2.out') 
    residual2_Ne6  = np.loadtxt('Newton/residual2(t)_Ne6_tol1e-2.out') 
    residual3_Ne6  = np.loadtxt('Newton/residual3(t)_Ne6_tol1e-2.out')  
    
    residual_directsim_Ne7 = np.loadtxt('Newton/residual_sim(t)_Ne7_tol1e-2.out')
    residual1_Ne7 = np.loadtxt('Newton/residual1(t)_Ne7_tol1e-2.out') 
    residual2_Ne7  = np.loadtxt('Newton/residual2(t)_Ne7_tol1e-2.out') 
    residual3_Ne7  = np.loadtxt('Newton/residual3(t)_Ne7_tol1e-2.out')  
#    
    Dtlist = [ 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1]
    
    
    ax = plt.subplot(111)  
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])   
   # Nlist = scipy.array([1e4,5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7]) 
       
   # fig = plt.figure()
   # ax = fig.add_subplot(1,1,1)
 #   ax = plt.subplot(111)  
  #  box = ax.get_position()
  #  ax.set_position([box.x0, box.y0, box.width * 0.82, box.height])   
          
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(Dtlist, residual1, "g--")
    plt.plot(Dtlist, residual2, "r--" )
    plt.plot(Dtlist, residual3, "y--" )
    plt.plot(Dtlist, residual_directsim, "c--" )
    
    plt.plot(Dtlist, residual1_Ne6, "g-.")
    plt.plot(Dtlist, residual2_Ne6, "r-." )
    plt.plot(Dtlist, residual3_Ne6, "y-." )
    plt.plot(Dtlist, residual_directsim_Ne6, "c-." )    
    
    plt.plot(Dtlist, residual1_Ne7, "g-", label=r'$k_{Newton} = 1$')
    plt.plot(Dtlist, residual2_Ne7, "r-", label = r'$k_{Newton} = 2$')
    plt.plot(Dtlist, residual3_Ne7, "y-",  label= r'$k_{Newton} = 3$')
    plt.plot(Dtlist, residual_directsim_Ne7, "c-", label = 'direct sim' )    
#    
#    
    plt.xlabel(r'$\Delta T$', fontsize = 14)
   # plt.xlabel(r'$N$', fontsize = 14)
    plt.ylabel(r'$||  \rho^{*}   - \frac{ \exp{\left[-2\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}  ||$'  , fontsize = 15)
    #plt.legend(prop={'size':6})
    
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2,3,4,5)
    #Create custom artists
    N5Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
    N6Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='-.')
    N7Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='-')
    
    ax.legend( numpoints = 1,  prop={'size':8} )
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[N5Artist,N6Artist, N7Artist],
          [label for i,label in enumerate(labels) if i in display]+[r'$N = 10^{5}$', r'$N = 10^{6} $',  r'$N = 10^{7} $'], bbox_to_anchor=(1.45
          , 1), prop={'size':9})    
    
    plt.savefig("plots/Newton_sde_res(Dt)_dt=1e-3_tol_1e-2.pdf")
    plt.show()
#    plt.savefig("plots/Newton_sde_res(it)_Dt=1e-1_tol_1e-3_alpha99_dxfine.pdf")   
#    plt.plot(resnorm_t10_tol2_N1000,  "g-", label=r'$N=10^3$')
#    plt.plot( resnorm_t10_tol2_N10000, "r-",label = r'$N=10^4$')
#    plt.plot( resnorm_t10_tol2_N100000,"y-" , label= r'$N=10^5$')
#    plt.plot(  resnorm_t10_tol2_N1000000, "c-" , label = r'$N=10^6$')    
#    plt.plot(  resnorm_t10_tol2_N10000000, "b-" , label = r'$N=10^7$')   

#    plt.plot( resnorm_t10_N1000,  "g--")
#    plt.plot( resnorm_t10_N10000, "r--")
#    plt.plot( resnorm_t10_N100000,"y--" )
#    plt.plot( resnorm_t10_N1000000, "c--"  )    
#    plt.plot(  resnorm_t10_N10000000, "b--" )   
#    
#  #  plt.xlabel(r'$\Delta T$', fontsize = 14)
#    plt.xlabel(r'$k$', fontsize = 14)
#    plt.ylabel(r'$||  \rho - \Phi^N_T(\rho)||$'  , fontsize = 15)
#
#  #  plt.title(r'$\Delta T = 0.01 , tol(-) = 10^{-2}, tol(--)=10^{-3} $')
#    
#    #Get artists and labels for legend and chose which ones to display
#    handles, labels = ax.get_legend_handles_labels()
#    display = (0,1,2,3,4)
#    #Create custom artists
#    tol3Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='--')
#    tol2Artist = plt.Line2D((0,1),(0,0), color='k', linestyle='-')
#    
#    ax.legend( numpoints = 1,  prop={'size':8} )
#    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[tol3Artist,tol2Artist],
#          [label for i,label in enumerate(labels) if i in display]+[r'$\epsilon_{GMRES} = 10^{-3}$', r'$\epsilon_{GMRES} = 10^{-2}$'], bbox_to_anchor=(1.39
#          , 1), prop={'size':10})
#  
#   # ax.legend(prop={'size':6})n
#   # plt.savefig("plots/Newton_sde_res(Dt)_N5e5dt=1e-3_tol_1e-3.pdf")
#    plt.savefig("plots/Newton_sde_res(k)_Dt_e-2_alpha99.pdf")   
#
#        


        

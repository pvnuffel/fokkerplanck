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


    res_e3= np.loadtxt('GMRES/res_Dt_e-2_builtin_N1000.out') 
    res_e4= np.loadtxt('GMRES/res_Dt_e-2_builtin_N10000.out') 
    res_e5= np.loadtxt('GMRES/res_Dt_e-2_builtin_N100000.out') 
    res_e6= np.loadtxt('GMRES/res_Dt_e-2_builtin_N1000000.out') 

    Nlist = scipy.array([10**3, 10**4, 10**5, 10**6])
       
    plot_res= plt.plot(res_e3, 'yo' , label= r'$N=10^3$', markersize =5)# % Nlist[0])  #gives the smae result, as expected
    plot_res= plt.plot(res_e4, 'bo' ,label= r'$N=10^4$', markersize =5) #)% Nlist[1]) 
    plot_res= plt.plot(res_e5, 'ro' ,label= r'$N=10^5$' , markersize =5)#% Nlist[2]) 
    plot_res= plt.plot(res_e6, 'go' , label= r'$N=10^6$', markersize =5) #% Nlist[3]) 
        
    plt.ylabel('GMRES residual' , fontsize = 11)
    plt.xlabel('iterations', fontsize = 11)
    plt.legend([plot_res], loc='best') 
    plt.legend(bbox_to_anchor=(1, 1), numpoints = 1 )    
    
     #   plt.xscale(log_flag)    res_e3= np.loadtxt('GMRES/res_N10000.out') 
    
    plt.yscale('log')
    plt.savefig("plots/GMRES/GMRES_N_Dt_1e-2_builtin.pdf")
       
    
    plt.show()
    
    
    
    
    spectrum_e3= np.loadtxt('GMRES/spectrum_Dt_e-2_builtin_N1000.out').view(complex) 
    spectrum_e4= np.loadtxt('GMRES/spectrum_Dt_e-2_builtin_N10000.out').view(complex) 
    spectrum_e5= np.loadtxt('GMRES/spectrum_Dt_e-2_builtin_N100000.out').view(complex)  
    spectrum_e6= np.loadtxt('GMRES/spectrum_Dt_e-2_builtin_N1000000.out').view(complex)  
    
    plot_spectrum= plt.plot(spectrum_e3.real, spectrum_e3.imag,  'yo' , label= r'$N=10^3$', markersize =4)
    plot_spectrum= plt.plot(spectrum_e4.real, spectrum_e4.imag,  'bo' , label= r'$N=10^4$', markersize =4)
    plot_spectrum= plt.plot(spectrum_e5.real, spectrum_e5.imag,  'ro' , label= r'$N=10^5$', markersize =4)
    plot_spectrum= plt.plot(spectrum_e6.real, spectrum_e6.imag,  'go' , label= r'$N=10^6$', markersize =4)
    
    plt.legend([plot_res], loc='best') 
    plt.legend(bbox_to_anchor=(0.3, 1), numpoints = 1 )  
    plt.savefig("plots/GMRES/Spectrum_Dt_1e-2_builtin.pdf")
#    plot_res= plt.plot(res_e4, 'bo' ,label= r'$N=10^4$', markersize =5) #)% Nlist[1]) 
#    plot_res= plt.plot(res_e5, 'ro' ,label= r'$N=10^5$' , markersize =5)#% Nlist[2]) 
#    plot_res= plt.plot(res_e6, 'go' , label= r'$N=10^6$', markersize =5) #% Nlist[3]) 
#    
    
            


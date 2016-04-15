"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np

from scipy.linalg import norm
from scipy import sqrt

#import matplotlib.pylab as plt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")

import pde
import histogram
import particles
import inv_transform

if __name__=="__main__":
    D = 1./2.
    Dt = 1e-2
    seed = 16
    xL = -1.7
    xR = 1.7
    
    a = 1
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([a,D,alpha,beta,zeta])
  
    #SDE
    t1 = time.time()
    
    print 'Start solving sde'
    dx = 1e-2
    dt = 1e-4 # ti
    N=1000
    Nlarge= N
    Nsmall =N
    
    
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    
    

    v=scipy.zeros_like(grid)
    for j in range(len(grid)): 
        v[j]= np.sin(j*2*np.pi/len(grid))

     
    
  #  sampler_param = inv_transform.Sampler.getDefaultParameters()
  #  sampler_param['seed']=0
  #  lifting =  inv_transform.Sampler(sampler_param)
    
    M=100 #number of monte carlo steps                

    restriction = histogram.Histogram()

    param=particles.Particles.getDefaultParameters()
    param['Nlarge']=N
    param['Nsmall']=N    
    param['Dt'] = Dt 
    param['dt']=dt
#    param['eps']=1e-5
   # fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param=param)       

    Nlist = scipy.array([1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000])
 #   Nlist = scipy.array([1000, 2000, 4000])  #, 800, 1600, 3200])
    eps_list_exponents=[4,5,6,7]
    eps_list= [1e-4, 1e-5, 1e-6, 1e-7]
   # Nlist = scipy.array([5,10, 200]) #,400])
  #  Nlist = scipy.array([4000])   

    nN= len(Nlist)
    
    
    for eps_i in range(len(eps_list)):
        
        param['eps']= eps_list[eps_i]
        print "----------------------------  \n",
        print "----------------------------  \n",
        print "running with epsilon ", eps_list[eps_i]
    
        rho_sq = scipy.zeros(nN)
        Jv_sq = scipy.zeros(nN)
        sq_E_rho=scipy.zeros(nN)
        sq_E_Jv=scipy.zeros(nN)
        E_rho= scipy.zeros((nN,len(rho)))
        E_Jv= scipy.zeros((nN,len(rho)))
    
        for m in range(1,M+1):    #M steps
            print "running with seed ", m 
            sampler_param = inv_transform.Sampler.getDefaultParameters()
            sampler_param['seed']=m
            lifting = inv_transform.Sampler(sampler_param)         
            #fp_sde=None
            x_prev_sim= None
            w_prev_sim = None
            rg_state = None
            
            for n in range(nN):  
                N = Nlist[n] #- np.sign(n)*Nlist[n-1]   #reduce the number of particle simulations by using the data from previous simulations
                if(len(grid)>N):
                    print "The number of particles is smaller than the number of bins! Please increase N or increase dx"
    
               # param['Nlarge']=Nlist[n]
                param['N']=N    
               # param['Nsmall']=Nlist[n]    
            
                print 'run simulation with N = ', N, 'particles'
                if(rg_state != None): 
                    lifting.set_rg_state(rg_state)
                fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param=param )   
             #   rg_state = lifting.get_rg_state()  #We remember the internal state of the random generator to get new random number for sampling the new particles
                
         #       print "calculate .u_Dt"
                rho_Dt = fp_sde.u_Dt 
                #rho_fine = fp_sde.make_rho_fine(rho_Dt) 
                rho_sq[n] = rho_sq[n] + (norm(rho_Dt))**2
                E_rho[n] = E_rho[n] + rho_Dt           #matrix
        
                print "calculate JVs"
                Jv =fp_sde.applyJacobian(v)                                                 
                Jv_sq[n] = Jv_sq[n] + (norm(Jv))**2
                E_Jv[n] = E_Jv[n] + Jv 
                
          #      x_prev_sim = fp_sde.x_Dt
              #  print 'len x-vector= ' , len( fp_sde.x_Dt) 
              #  w_prev_sim = fp_sde.w_prev
                #print w_prev_sim
                                 
        for n in range(len(Nlist)): 
            rho_sq[n]=  rho_sq[n]/M
            E_rho[n] = E_rho[n]/M
            sq_E_rho[n]= (norm(E_rho[n]))**2
            
            Jv_sq[n]=  Jv_sq[n]/M
            E_Jv[n] = E_Jv[n]/M
            sq_E_Jv[n]= (norm(E_Jv[n]))**2
            
        
        np.savetxt('nw-Jv_sq_eps-%d.out' %eps_list_exponents[eps_i], Jv_sq)

        np.savetxt('nw-sq_E_Jv_eps-%d.out' %eps_list_exponents[eps_i], sq_E_Jv)
    

    print "End of simulation"
    now = time.time()
  #  print "Simulation time for solving pde: " , t1-t0, " seconds"
    print "Simulation time for solving sde: " , now-t1, " seconds"
#    print "Total simulation time " , now-t0, " seconds"
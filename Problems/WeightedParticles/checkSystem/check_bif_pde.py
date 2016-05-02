# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:05:37 2015

@author: pieter
"""
import time 
import tables
import scipy
import sys
import numpy as np
from pylab import *
import os


from scipy.linalg import norm
from scipy import sqrt

import matplotlib.pylab as plt

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../..")

import pde

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability
import Continuer.Continuer as Continuer 





if __name__=="__main__":
   # D = 1./2.
    mu = 1.0
    sigma=2.0
    D=0.5*sigma**2
  #  D = 2.0
   # sigma = 2*np.sqrt(D)
    Dt = 1e-2
  #  Dtlist = [1e-5,2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3 , 1e-2, 2e-2, 5e-2 , 1e-1, 2e-1, 5e-1, 1, 2, 5, 10]
    Dtlist = [Dt]
 #   Dtlist = [2e-3]
    # discretizatio    seed = 15n parameters
   #h=2e-2 # kde mesh width 
    dx = 1e-2
#    dxlist= [1e-2,  5e-3, 1e-3, 5e-4] 
  #  dxlist= [1e-2, 1e-3]
    #dt = min(dxlist)**2/(2*D)
    dt=1e-5
    r= (2.*D*dt)/(dx)**2
    print'r = ' ,  r

   # r= (2.*D*dt)/(min(dxlist))**2
    if r>1: 
        print 'Stability condition not fulfilled, because r=',r #(see https://en.wikipedia.org/wiki/FTCS_scheme)
        sys.exit("Ending Script. Make sure that r<1") 
    xL = -1.7
    xR = 1.7
    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([mu, sigma])
    param=pde.FokkerPlanck.getDefaultParameters()
    param['Dt'] = Dt
    param['dt']=dt         
  #  norm_Jv_fp= scipy.zeros(len(dxlist))
 #   norm_rho_fp= scipy.zeros(len(dxlist))
   # rho_Dt_fp[i]= np.ndarray(shape=) #[scipy.zeros((len
    

    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho
    
    residual = scipy.zeros(len(Dtlist))
    residual2 = scipy.zeros(len(Dtlist))
    residual_directsim = scipy.zeros(len(Dtlist)) 
    t0=time.time()
    param['eps']=1e-5   #strongly influences number of gmres-iterations 
    
    for dti in range(0, len(Dtlist)):
        Dt = Dtlist[dti]
        param['Dt']=Dt       
        
        print 'Total Simulation Time Dt = ' , Dt
        fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
        rho_Dt=fp_pde.u_Dt
    
    #    v=scipy.zeros_like(grid)
    # #   v[200]=-1
    # #   v[300]=1
    #    for j in range(len(grid)): 
    #        v[j]= np.sin(j*2*np.pi/len(grid))
    
        #  plt.show()    
        # CREATING LINEAR SOLVER
        gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
        gmres_param['tol']=1e-5
        gmres_param['print']='short'
        gmres_param['builtin']=True
        linsolv = GMRES.GMRESLinearSolver(gmres_param)
        linsolv2 = GMRES.GMRESLinearSolver(gmres_param)
    
        # CREATING NEWTON SOLVER
        newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
        newt_param['rel_tol']=1e-5
        newt_param['abs_tol']=1e-6
        newt_param['print']='short'
        newt_param['max_iter']=5
        newt_param['damping']=1
        nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
     #   nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)
    
        # CREATING POINT
        psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
        psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolver()
    
        # POINT PARAMETERS
        point_param = Point.Point.getDefaultParameters()
      #  point_param['free']= 
     #   point_param['artificial']=[2]
     #   point_param['artificial_condition']=[probability.condition]
        
     #   pprev = p
        points = []
        p = Point.Point(fp_pde,nsolv,None, point_param)
        print "Solve for fixpoint"
        p.correct()
        print "Found! Add to list"
        points.append(p)        
               
       # spectrum = p.mspectrum()
                
        newton_states = nsolv.newton_states
        nstates =  newton_states.reshape(nsolv.nb_newt+1, len(grid))
        
        norm_c = sum( np.exp( (-grid**4 + grid**2)*mu/D))*dx
        rho_ss = np.exp( (-grid**4 + grid**2 )*mu/D) /norm_c
        label2= r'$\frac{ \exp{\left[-\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}$'
        residual_directsim[dti] = norm( rho_ss- rho_Dt)
        residual[dti] = norm( nstates[-1] - rho_ss)
       # residual2[dti] = norm( nstates[-2]-rho_ss)
        
        
#        lambd2 = scipy.array([a,0.4,alpha])
#        fp_pde2 = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd2,param)
#        p2 = Point.Point(fp_pde2,nsolv,None, point_param)
#        p2.correct()
#        print "Found! Add to list"
#        points.append(p2)     
#        
        continuer_param= Continuer.Continuer.getDefaultParameters()
        continuer_param['plotx']['func'] = 1
        continuer_param['growth_factor']=1
        continuer_param['cont_step']=-0.1
        branch_pde = Continuer.Continuer(points, continuer_param)
   #     branch.bcontinue(4)
        print '#Start Continuation'
        cont_steps=10
        succ = branch_pde.bcontinue_natural(cont_steps)
        x_pde, y_pde = branch_pde.getPlot()
        N_iter_pde = branch_pde.getIterations()
        print succ, 'succesfull continuation steps', cont_steps -succ, 'failed ones'
        
                    
        newton_res_norms = nsolv.newton_res_norm
     #   np.savetxt('Newton/pde_res_Dt1e-3.out', newton_res_norms)
     #   np.savetxt('GMRES/pde_Dt1e-3_spectrum_Dt', spectrum[0].view(float)) 
       # plot_spectrum = plt.plot(spectrum[0].real, spectrum[0].imag,  'go' , markersize =4)    
#        plt.show()

        for i in range(nstates.shape[0]):
            plot_rho = plt.plot(grid, rho_Dt, 'gray')
            plot_rho = plt.plot(grid, nstates[i], "g-",  linewidth=2 )
            plot_rho = plt.plot(grid, rho_ss,  "r--",  linewidth=2, label=label2)
         #   plot_rho = plt.plot(grid, rho_Dt[i], "c--", linewidth=2)
            plot_rho = plt.plot(grid, nstates[i]-rho_ss, "y--", linewidth=2)
            plt.show()
            

    
  
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(Dtlist, residual, "g-")
  #  plt.plot(Dtlist, residual2, "r--" )
    plt.plot(Dtlist, residual_directsim, "c--" )
    plt.xlabel(r'$\Delta T$', fontsize = 14)
    plt.ylabel(r'$||  \rho^{*}   - \frac{ \exp{\left[-\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}  ||$'  , fontsize = 15)
   # plt.savefig("plots/Newton_pde_res(Dt)_dt=1e-5.pdf")
              
    
    
#    sigma =2
#    lambd = scipy.array([a,sigma,alpha,beta,zeta])
#    fp_pde2 = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)
#    p2 = Point.Point(fp_pde2,nsolv2,None,point_param)
#    p2.correct()
#    points.append(p2)            
#      

#           
#
#    
    
        #rho_Dt_fp = fp_pde.u_Dt
        #plt.plot( grid, rho_Dt_fp)
     #   plt.plot( rho_Dt_fp)  
     #   norm_rho_fp[i] = norm(rho_Dt_fp)/sqrt(len(rho_Dt_fp))
       # norm_Jv_fp[i] = norm(Jv_pde)/sqrt(len(Jv_pde))
        
        
    print "End of simulation"
    now = time.time()
    print "Simulation time for calculating fixpoint: " , now-t0, " seconds"
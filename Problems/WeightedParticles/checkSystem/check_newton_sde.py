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

import histogram
import particles
import inv_transform_sampling as inv_transform

#import matplotlib.pylab as plt

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.Richardson as Rich
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability
#import Continuer.Continuer as Continuer 
import Continuer.Continuer as Continuer 


if __name__=="__main__":
    print 'Start solving sde'
 #   sigma=2
    mu = 1.0
  #  D=0.5*sigma**2
    D = 2.0
    sigma=1.0
   # sigma = 2*np.sqrt(D)
    Dt = 1e-2  #Dt influences number of GMRES-iterations strongly
    Dtlist = [ Dt]
    #Dtlist = [ 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1]
    #Dtlist = [1e-5, 5e-5, 1e-4,  5e-4, 1e-3,  5e-3 , 1e-2, 5e-2 , 1e-1, 5e-1, 1]
    seed = 1
    # discretization parameters
   #h=2e-2 # kde mesh width 
   # dt=1e-6
    xL = -1.7
    xR = 1.7

    zeta = 0.
    alpha=1
    beta = -1
    lambd = scipy.array([mu, sigma])
       
    #SDE
    t1 = time.time()
    
    print 'Start solving sde'
    dx = 1e-2
    dt = 1e-3 # ti
#    N=5e5
#    Nlarge= N 
#    Nsmall =N
    
    grid = scipy.arange(xL+dx/2.  ,xR,dx)
    rho = scipy.ones_like(grid)
    rho = rho/(sum(rho)*dx)
#    v=scipy.zeros_like(grid)
#    for j in range(len(grid)): 
#        v[j]= np.sin(j*2*np.pi/len(grid))
    
  #  sampler_param = inv_transform.Sampler.getDefaultParameters()
  #  sampler_param['seed']=0
  #  lifting =  inv_transform.Sampler(sampler_param)
    
    h=2e-2 # kde mesh width     
    M=1 #number of monte carlo steps         
    y_mean = scipy.zeros(M)
    y_sde = scipy.zeros(M)

    branch_list = [] #scipy.array([M])              
    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=h
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
   # param['Nlarge']=N
   # param['Nsmall']=N    
    param['Dt'] = Dt 
    param['dt']=dt
#    param['eps']=1e-5
   # fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, param=param)       

 #   Nlist = scipy.array([1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000])
    Nlist = scipy.array([ 1e7]) 
  #  Nlist = scipy.array([1e5]) #, 800, 1600, 3200])
    eps_list_exponents=[5]
    eps_list= [1e-5]
    param['eps']= eps_list[-1]
   # Nlist = scipy.array([5,10, 200]) #,400])-
  #  Nlist = scipy.array([4000])   
    
    nN= len(Nlist)
    
    rho_sq = scipy.zeros(nN)
    Jv_sq = scipy.zeros(nN)
    sq_E_rho=scipy.zeros(nN)
    sq_E_Jv=scipy.zeros(nN)
    E_rho= scipy.zeros((nN,len(rho)))
    E_Jv= scipy.zeros((nN,len(rho)))
    
   # = scipy.zeros((nN,len(rho)))
    
    
        
#    residual1 = scipy.zeros(len(Dtlist))
#    residual2 = scipy.zeros(len(Dtlist))
#    residual3 = scipy.zeros(len(Dtlist))
#    residual_directsim = scipy.zeros(len(Dtlist)) 
    save_flag= False

#            
#    residual1 = scipy.zeros(len(Dtlist))
#    residual2 = scipy.zeros(len(Dtlist))
    residual = scipy.zeros(11)

    residual_directsim = scipy.zeros(nN) 
    newton_res_norms = scipy.zeros(nN)
    

    for Dti in range(0, len(Dtlist)):
        Dt = Dtlist[Dti]
        param['Dt']=Dt       
        print 'Total Simulation Time Dt = ' , Dt    
          
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
               # N = Nlist[n] - np.sign(n)*Nlist[n-1]   #reduce the number of particle simulations by using the data from previous simulations
                N = Nlist[n]
                if(len(grid)>N):
                    print "The number of particles is smaller than the number of bins! Please increase N or increase dx"
    
    #            param['Nlarge']=Nlist[n]
               
                param['Nsmall']=N    
                param['Nlarge']=N    
            
                print 'run simulation with N = ', N, ' = ', Nlist[n], ' - ' , np.sign(n)*Nlist[n-1]  , 'particles'
                if(rg_state != None): 
                    lifting.set_rg_state(rg_state)
                fp_sde = particles.Particles(lifting,restriction,rho,grid, lambd, x_prev_sim, w_prev_sim, param=param )   
                rg_state = lifting.get_rg_state()  #We remember the internal state of the random generator to get new random number for sampling the new particles
                
                print "calculate .u_Dt"
             #   rho_Dt = fp_sde.u_Dt 
                t0=time.time()
                print "Simulation time for solving sde: " , t0-t1, " seconds"
                #rho_fine = fp_sde.make_rho_fine(rho_Dt) 
             #   rho_sq[n] = rho_sq[n] + (norm(rho_Dt))**2
                    #matrix
        
                print "calculate JVs"
              #  Jv =fp_sde.applyJacobian(v)                                                 
               # Jv_sq[n] = Jv_sq[n] + (norm(Jv))**2
                #E_Jv[n] = E_Jv[n] + Jv 
                
                            #LINEAR SOLVER
                gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
                gmres_param['tol']=1e-1   #zet terup 1e-4
                gmres_param['print']='none'
                gmres_param['builtin']=False
                
                linsolv = GMRES.GMRESLinearSolver(gmres_param)
              #  richardson_param = Rich.Richardson.getDefaultParameters()
              #  linsolv = Rich.Richardson(richardson_param)
              #  linsolv2 = GMRES.GMRESLinearSolver(gmres_param)
        
                # CREATING NEWTON SOLVER
                newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
                newt_param['rel_tol']=1e-5
                newt_param['abs_tol']=1e-10
                newt_param['print']='short'
                newt_param['max_iter']=8
                newt_param['damping']=1
                nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
            #    nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)
                # CREATING POINT
        #        psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
             #   psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolv@er()
        
                # POINT PARAMETERS
                point_param = Point.Point.getDefaultParameters()
             #   point_param['artificial']=[4]
             #   point_param['artificial_condition']=[probability.condition]
             #   pprev = p
                p = Point.Point(fp_sde,nsolv,None, point_param)
                points = []
                print "Solve for fixpoint"
                p.correct()
              #  print "Found! Add to list"
               # points.append(p)       
         
               
                        

              #  x_prev_sim = fp_sde.x_Dt
              #  print 'len x-vector= ' , len( fp_sde.x_Dt) 
             #   w_prev_sim = fp_sde.w_prev
                #print w_prev_sim
                

#                p = Point.Point(fp_sde,nsolv,None, point_param)   #if using  (psolver-im)-> ATTENTION: RUNNING COMPUTE JACOBIAN ON TIMESTEPPER SYSTEM
#                p.correct()stabiel systeem
#                
              #  spectrum = p.mspectrum()
                
           #     rho_fix =p.u 
           #     E_rho[n] = E_rho[n] +   rho_fix      
             #   rho_sq[n] = rho_sq[n] +   norm(rho_fix)**2
               # kappa= abs(scipy.amax(spectrum[0]))/abs(scipy.amin(spectrum[0]))
                # points.append(p)
                newton_states = nsolv.newton_states
                #newton_res_norms = nsolv.newton_res_norm
             #   np.savetxt('Newton/new_method_resnorm_Dte-2_tole-4_Ne6m%d' %m, newton_res_norms)
                #np.savetxt('Newton/21-05_resnorm_Dte-2_tole-7_N%d' %N, newton_res_norms)
                States =  newton_states.reshape(nsolv.nb_newt +1, len(grid))  
                norm_c = sum( np.exp( 2*(-grid**4 + grid**2)*mu/sigma**2))*dx
                rho_ss = np.exp( 2*(-grid**4 + grid**2 )*mu/sigma**2) /norm_c 
                #label2= r'$\frac{ \exp{\left[-2\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}$'         
            #    np.savetxt('GMRES/res_Dt1e-1_dx1e-1_N%d.out' %N, linsolv.resid)
            #    print linsolv.resid
            #    np.savetxt('GMRES/spectrum_Dt1e-1_dx1e-1_N%d.out' %N, spectrum[0].view(float))
              #  np.savetxt('Newton/states.out', nsolv.newton_states)
        
             #   residual_directsim[Dti] = norm( rho_Dt -  rho_ss)/sqrt(len(rho_Dt- rho_ss))
             #  residual2[Dti] = norm( States[2] - rho_ss)/sqrt(len(States[2] - rho_ss))
             #   residual1[Dti] = norm(States[1]-rho_ss)/sqrt(len(States[1] - rho_ss))
             #   residual3[Dti] = norm(States[3]-rho_ss)/sqrt(len(States[3] - rho_ss))
        
        #        residual_directsim[n] = norm( rho_Dt -  rho_ss)/sqrt(len(rho_Dt- rho_ss))
                for k in range(len(States)):
                    residual[k] = norm(States[k]-rho_ss)/sqrt(len(States[k] - rho_ss))
                    
           #     np.savetxt('Newton/25-05_nstates_t100/25_05_residual__tol-4_t_100_N%d' %N, residual)
                np.savetxt('Newton/Variable_GMREStol/residual_8it_N1e7',  residual)
                np.savetxt('Newton/Variable_GMREStol/Newton_states_8it_N1e7', States)
           #            np.savetxt('Newton/26-05_Dependency_of_GMREStol_N1e6/Newton_states_tol-6', States)

#                for i in range(States.shape[0]):
#                    print i
#                 #   plot_rho = plt.plot(grid, States[i], "g-",  linewidth=2 )
#                    plot_rho = plt.plot(grid, rho_ss,  "r--",  linewidth=2)
#                    if (i>1):
#                        rho_mean = rho_mean + States[i]
#                        plot_rho_mean = plt.plot(grid, rho_mean/(i-1),  "b-",  linewidth=1)
#                    plot_rho = plt.plot(grid, p.u,  "g-",  linewidth=2)
#                 #   plot_rho = plt.plot(grid, rho_Dt, "c--", linewidth=2)
#              #      label = 'Newton/Newtonstep%d.pdf' %i
#                    #plt.savefig(label) 
#                    plt.show()
#    
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.plot(Dtlist, residual1, "g--")
#    plt.plot(Dtlist, residual2, "r--" )
#    plt.plot(Dtlist, residual3, "y--" )
#    plt.plot(Dtlist, residual_directsim, "c--" )
#    plt.xlabel(r'$\Delta T$', fontsize = 14)
#    plt.ylabel(r'$||  \rho^{*}   - \frac{ \exp{\left[-2\frac{V(x)}{\sigma^2}\right]}}{\mathcal{N}}  ||$'  , fontsize = 15)
#    plt.savefig("plots/Newton_sde_res(Dt)_dt=1e-5.pdf")
              
#   np.savetxt('Newton/BuiltinGMRES_tole-2_N1e5/residual_sim(N)_Dte-1.out' , residual_directsim)
#    np.savetxt('Newton/residual_sim(t)_Ne7_tol1e-2.out' , residual_directsim)
#    np.savetxt('Newton/residual1(t)_Ne7_tol1e-2.out' , residual1)
#    np.savetxt('Newton/residual2(t)_Ne7_tol1e-2.out' , residual2)
#    np.savetxt('Newton/residual3(t)_Ne7_tol1e-2.out' , residual3)
              
#    for n in range(len(Nlist)): 
#        rho_sq[n]=  rho_sq[n]/M
#        E_rho[n] = E_rho[n]/M
          #  sq_E_rho[n]= (norm(E_rho[n]))**2
            
          #  Jv_sq[n]=  Jv_sq[n]/M
          #  sq_E_Jv[n]= (norm(E_Jv[n]))**2

   # np.savetxt('Newton/rho_sq_tol_e-5.out' , rho_sq)
   # np.savetxt('Newton/E_rho_tol_e-5.out' , E_rho )


    
##       
    print "End of simulation"
    now = time.time()
    print "Simulation time for calculating fixpoint: " , now-t0, " seconds"
    print "Simulation time for solving sde: " , now-t1, " seconds"
#    print "Total simulation time " , now-t0, " seconds"
import sys

sys.path.append('../')
sys.path.append('../../')
sys.path.append('../../../')
sys.path.append('../system/')
sys.path.append('../lifting/')
sys.path.append('../restriction/')
sys.path.append('../modules/')

import scipy
import scipy.linalg

import inv_transform
#import kde
import histogram
import particles
import pde
import precond

# Fokker Planck for SDE
a = 1  #advection coefficient
D = 0.1  #difussion coefficient
lambd = scipy.array([a,D])
Dt = 0.01 
N = 100000
eps = 1e-7

# discretization parameters
h=2e-2 # kde mesh width 
dx = 0.05
dt = 1e-3 
print("check stability condition")
print ("nu=", dx**2/(2.*D*dt))  #stability condition: nu>1  (see https://en.wikipedia.org/wiki/FTCS_scheme)
xL = -2.
xR = 2.
grid = scipy.arange(xL+dx/2.,xR,dx)
rho = scipy.ones_like(grid)/(xR-xL)
print ("rho: ", rho)


# Fokker Planck for SDE
sampler_param = inv_transform.Sampler.getDefaultParameters()
sampler_param['seed']=0
lifting = inv_transform.Sampler(sampler_param)

# param_kde = kde.KDE.getDefaultParameters()
# param_kde['h']=h
# restriction = kde.KDE(param_kde)

restriction = histogram.Histogram()

param=particles.Particles.getDefaultParameters()
param['eps']=eps
param['Dt'] = Dt
param['N']=100000    
precond_param=precond.Precond.getDefaultParameters()
precond_param['Dstar']=0.
precond_param['sigma']=scipy.zeros_like(grid)
precond_param['kappa']=particles.doublewell
param['precond']=precond.Precond(precond_param)

fp_sde = particles.Particles(lifting,restriction,rho,grid,lambd,param)

v = scipy.zeros_like(grid)
v[10]=1.
v[20]=-1.

# testing the run
print (" run  ----")
rho_Dt = fp_sde.integrate(rho,lambd)
print (" perturbed run ---")
rho_Dt_pert = fp_sde.integrate(rho+eps*v,lambd)
print ("Jacobian --- ")
Jv_sde = fp_sde.applyJacobian(v)

# # testing the sampler algorithm :
# 
# rho_Dt = scipy.zeros_like(grid)
# rho_Dt[5:-5] = 1.
# rho_Dt = rho_Dt/sum(rho_Dt)*(grid[1]-grid[0])
# xx = sampler.sample_density(rho_Dt,grid,N)
# rho_Dt_resampled = fp_sde.restrict_histogram(xx)

# Corresponding PDE

param=pde.FokkerPlanck.getDefaultParameters()
param['D']=D
param['Dt'] = Dt
param['eps_jac']=eps

fp_pde = pde.FokkerPlanck(rho,grid,pde.doublewell,lambd,param)

Jv_pde = fp_pde.applyJacobian(v)

# testing the preconditioner

#A = scipy.identity(len(grid))+Dt*fp_sde.computeJacobian()
#J = fp_sde.computeJacobian()
#pJ = scipy.linalg.solve(A,v)

print (sum(v), sum(Jv_sde), sum(Jv_pde) )

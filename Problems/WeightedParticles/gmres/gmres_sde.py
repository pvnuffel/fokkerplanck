import tables
import scipy
import sys

sys.path.append("../system/")
sys.path.append("../lifting/")
sys.path.append("../restriction/")
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../..")

import Point.Point as Point
import Solver.NewtonSolver as NewtonSolver
import Solver.GMRESLinearSolver as GMRES
import Solver.ImplicitEulerDirectLinearSolver as ImDirect
import Utils.conditions.probability as probability

#import inv_transform_sampling as inv_transform
#import inv_transform_determ as inv_transform
import inv_transform_sampling as inv_transform
import histogram
import particles
import pde
import precond

if __name__=="__main__":
     
    D = None
    Dt = None
    N = None
    # deze manier van input parameters verwerken staat in Langtangen
    while len(sys.argv) > 1 :
        option = sys.argv[1];               del sys.argv[1]
        if option == "-Dt" :
            Dt = float(sys.argv[1]);         del sys.argv[1]
            print "#Dt: ", Dt, " (from command line)"
        elif option == "-N" :
            N = int(sys.argv[1]);         del sys.argv[1]
            print "#N: ", N, " (from command line)"
            
    if D == None:
        D = 1./2.
    if Dt == None:
        Dt = 1e-3
    if N == None:
        Nlarge = 100000
        Nsmall = 1000
        
    seed = 16
    # discretization parameters
    h=2e-2 # kde mesh width 
    dx = 0.05
    dt = 1e-3
    print "nu : ", dx**2/2./D/dt
    xL = -1.7
    xR = 1.7
    grid = scipy.arange(xL+dx/2.,xR,dx)
    rho = scipy.ones_like(grid)+\
        0.1*scipy.random.uniform(low=-0.5,high=0.5,size=scipy.shape(grid)[0])
    rho = rho/(sum(rho)*dx)
    print "rho: ", rho

    # Fokker Planck for SDE
    a = 1
    zeta = 0.
    alpha = 1.
    beta = -1

    p = None
    fp_sde = None

    # Fokker Planck for SDE
    lambd = scipy.array([a,D,alpha,beta,zeta])

    sampler_param = inv_transform.Sampler.getDefaultParameters()
    sampler_param['seed']=seed
    lifting = inv_transform.Sampler(sampler_param)
   
    param_histogram = histogram.Histogram.getDefaultParameters()
    param_histogram['h']=h
    restriction = histogram.Histogram(param_histogram)

    param=particles.Particles.getDefaultParameters()
    param['Dt'] = Dt
    param['Nlarge']=Nlarge   
    param['Nsmall']=Nsmall   
    param['dt']=dt
#    precond_param=precond.Precond.getDefaultParameters()
#    precond_param['Dstar']=1./2.
#    precond_param['sigma']=scipy.zeros_like(grid)
#    precond_param['kappa']=scipy.zeros_like(grid)
#    param['precond']=precond.Precond(precond_param)

    param_histogram = histogram.Histogram.getDefaultParameters()
    restriction = histogram.Histogram(param_histogram)
    
    fprev = fp_sde
    x_prev_sim= None
    w_prev_sim = None
#    fp_sde = particles.Particles(lifting,restriction,rho,grid,\
#        lambd,control = fp_sde, param=param)
    fp_sde= particles.Particles(lifting,restriction,rho,grid, lambd, x_prev_sim, w_prev_sim, param=param )   
    print fp_sde
    # print fp_sde2,fp_sde2.control

    # CREATING LINEAR SOLVER
    gmres_param = GMRES.GMRESLinearSolver.getDefaultParameters()
    gmres_param['tol']=1e-8
    gmres_param['print']='short'
    gmres_param['builtin']=True
    linsolv = GMRES.GMRESLinearSolver(gmres_param)
    # linsolv2 = GMRES.GMRESLinearSolver(gmres_param)

    # CREATING NEWTON SOLVER
    newt_param = NewtonSolver.NewtonSolver.getDefaultParameters()
    newt_param['rel_tol']=1e-7
    newt_param['abs_tol']=1e-7
    newt_param['print']='short'
    newt_param['max_iter']=1
    nsolv = NewtonSolver.NewtonSolver(linsolv,newt_param)
    # nsolv2 = NewtonSolver.NewtonSolver(linsolv2,newt_param)

    # CREATING POINT
    psolver_im = ImDirect.ImplicitEulerDirectLinearSolver()
    # psolver_im2 = ImDirect.ImplicitEulerDirectLinearSolver()

    # POINT PARAMETERS
    point_param = Point.Point.getDefaultParameters()
    point_param['artificial']=[4]
    point_param['artificial_condition']=[probability.condition]
    
    pprev = p
    p = Point.Point(fp_sde,nsolv,None,point_param)
    p.correct()
    # p2 = Point.Point(fp_sde2,nsolv2,psolver_im2,point_param)

#    p.correct()

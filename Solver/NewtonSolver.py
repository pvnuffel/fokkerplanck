import scipy
import Solver
import LinearSolver
from scipy import zeros,dot,r_


import matplotlib.pylab as plt

class NewtonSolver(Solver.Solver):

    def __init__(self,linear_solver,parameters=None):
        """ 
        input: 
        =====
            linear_solver (LinearSolver) 
                contains the linear solver that will be used in 
                each Newton iteration
            parameters (dict) 
                look at the docstring of getDefaultParameters() 
                to find out which fields there are
        behaviour:
        =========
            This class implements a Newton solver that stops when the 
            maximum number of iterations has been reached, OR the 
            relative OR absolute tolerance have been reached.
        """
        Solver.Solver.__init__(self,parameters)
        if isinstance(linear_solver,LinearSolver.LinearSolver):
            self.linsolv=linear_solver
        else:
            raise TypeError, "input argument " + linear_solver \
                + " should be a linear solver"
        self.nb_newt = 0 
        self.newton_residual = zeros((0,))
        self.newton_res_norm = zeros(0,)
       # self.newton_states = zeros((  param['max_iter'],len( self.point.getCurrentGuess()) ))
        self.newton_states = zeros((  0 ))
      #  self.iterations  = 0

    def getDefaultParameters():
        """  
        Returns a dictionary of default parameters 
            'maxIter'       : 10
            'relTol'        : 1e-6
            'absTol'        : 1e-8
            'damping'       : 1 (means no damping)
            'print'         : 'short'  (one line per iteration)
                other possibilities
                'long' (residual in each iteration)
                'none' (prints nothing)
        """
        param = {}
        param['max_iter']=10
        param['rel_tol']=1e-6
        param['abs_tol']=1e-8
        param['print']='short'
        param['damping']=1
        param['stop']="res"
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)
        
    def solve(self):
        """
        Perform Newton iterations, starting from x0.  
        The method returns status
        Status means:
            0: method has converged
            1: method has used maximum number of iterations 
               without convergence
                """
        iter = 0
        max_iter = self.param['max_iter']
        abs_tol = self.param['abs_tol']
        rel_tol = self.param['rel_tol']
        print_it = self.param['print']
        alpha = self.param['damping']
        stop = self.param['stop']
        x = self.point.getCurrentGuess()
        self.newton_states = x
        res = self.point.getResidual()
        res_norm = scipy.linalg.norm(res)/len(res)
        self.newton_res_norm  =res_norm
        print 'res norm init ', self.newton_res_norm
        grid = self.point.system.grid
        grid_dx = grid[-1] - grid[-2]
        if stop == "step":
            stop_norm = 0
        elif stop == "res":
            stop_norm = res_norm
            
        #while (iter < max_iter and (stop_norm > abs_tol or iter ==0)): geprobeerd om iedentieke toestanden te vermijden in opl van pde, maar als het residu te klein is ivm de tolerantie, dan heeft een newtonstap geen zin want dan wordt dx toch 0
        while (iter < max_iter and (stop_norm > abs_tol)):
            iter += 1
            self.nb_newt += 1
            dx, status = self.linsolv.solve(-res)        # res = rho - rho.Dt
            if print_it == 'long':  
                plt.plot(res)                  
                plt.title(r'Res before newton step = $ \rho - \phi(\rho) $ ')
                plt.show()
                plt.plot(dx)
                plt.title(r'$dx$ ')
                plt.show()
            
            step_norm = scipy.linalg.norm(dx)
            step_avg = scipy.sum(dx)
            # Handling of errors and printing
            if not status == 0:
                print "An error has occurred in solving " +\
                    "the linear system"
            # print the iterations
            if not print_it=='none':
                print "# Newton_iter: ", iter, " | norm res  :", \
                    res_norm, "norm step : ", step_norm, "avg step: ", step_avg
#            if print_it=='long':
#                print "Newton Residual"
#                for i in range(len(res)):
#                    print i, res[i]
#                print "---------------------"
#                print "Newton Current Guess (before this iteration)"
#                for i in range(len(x)):
#                    print i, x[i]
#                print "---------------------"
#                print " Newton Step"
#                for i in range(len(x)):
#                    print i, dx[i]
#                print "---------------------"
#                print "Sum : ", scipy.sum(dx)
#                print "--------"
#                print "Newton Current Guess (after this iteration)"
#                for i in range(len(x)):
#                    print i, x[i]+alpha*dx[i]
#                print "---------------------"
    
           # print 'check sum=1? norm before = ', sum(x)*grid_dx
            x=x+alpha*dx
         #   print 'check sum(dx)=0? norm dx = ', sum(dx)*grid_dx
            #x = x/sum(x*grid_dx)
          #  print 'check sum 1? norm after = ', sum(x)*grid_dx
            self.newton_states = r_[self.newton_states,x]
            self.point.setCurrentGuess(x)             # =  self.point.setState() +  self.point.updateState()
            res = self.point.getResidual()
            res_norm = scipy.linalg.norm(res)/len(res)
            print 'resNorm ', res_norm
            self.newton_res_norm = r_[self.newton_res_norm,res_norm]
            print self.newton_res_norm
          #  self.point.getState()
            if stop == "step":
                stop_norm = scipy.linalg.norm(dx)
            elif stop == "res":
                stop_norm = res_norm
        if not print_it == 'none':
            print "Newton: final residual : ", res_norm
        self.iterations = iter
        self.point.iterations=iter
        if iter == max_iter:
            status = 1  
            print 'Failed Newton iteration'
        else:
            status = 0
        return status

    def setPoint(self,point):
        self.linsolv.setPoint(point)
        Solver.Solver.setPoint(self,point)
    
    def new_build(self):
        self.linsolv.new_build()       

import scipy

class Sampler(object):    
    def __init__(self,param=None):
        if param == None :
            param = Sampler.getDefaultParameters()
        self.param = param
        self.rand = scipy.random.mtrand.RandomState()
        self.seed(param['seed'])
    
    def getDefaultParameters():
        param = {}
        param['seed']=0
        param['nb_grid']=50
        return param
    getDefaultParameters=staticmethod(getDefaultParameters)
    
    def lift(self,rho,grid,N):
        print "HIER !"
        K = self.param['nb_grid']
        M = N/K
        cum,edges = self.rho2cum(rho,grid)
        print K
        y = 1./K*scipy.arange(0.5,K,1.)
        x = self.inv_cum(y,cum,edges)
        print y
        print x
        xx = scipy.repeat(x,M)
        print K,M
        return xx
        
    def rho2cum(self,rho,grid):
        ngrid = len(grid)
        cum = scipy.zeros((ngrid+1,))
        edges = scipy.zeros_like(cum)
        # grid values are box centers
        edges[0]=3./2.*grid[0]-1./2.*grid[1]
        edges[1:ngrid]=(grid[:ngrid-1]+grid[1:ngrid])/2.
        edges[-1]=3./2.*grid[-1]-1./2.*grid[-2]
        
        cum[0]=0.
        for n in range(ngrid):
            if rho[n]<0 :
                print "encountered a negative density : ", rho[n]
            cum[n+1]=cum[n]+rho[n]*(edges[n+1]-edges[n])
        print "estimating cumulative density : last value ", cum[-1], " (should be 1.)"
        if scipy.absolute(cum[-1]-1)>1e-7:
            print "cumulative: ", cum
        return cum,edges
    
    def inv_cum(self,y,cum,edges):
        x = scipy.zeros_like(y)
        for i in range(len(y)):
            x[i]=self.inv_cum_single_particle(y[i],cum,edges)
        return x
    
    def inv_cum_single_particle(self,yn,cum,edges):
        i = 0
        while (yn>cum[i]):
            i+=1
        m = (cum[i]-cum[i-1])/(edges[i]-edges[i-1])
        xn = (yn-cum[i-1])/m+edges[i-1]
        return xn
    
    def seed(self,s):
        print "seeding the density sampler with seed : ", s
        self.cur_seed = s
        self.rand.seed(self.cur_seed)
    
    def getBrownianIncrement(self,N):
        K = self.param['nb_grid']
        M = N/K
        dW = scipy.zeros((1,M))
        dW[0,:]=self.rand.normal(loc=0.,scale=1.,size=M)
        dWW = scipy.repeat(dW,K,axis=0)
        return dWW.flatten()
    


# Here we put a number of tests :
if __name__ == "__main__":
    
    # construction of initial grid and density
    xL = 0.
    xR = 20.
    Dx = 0.2
    grid = scipy.arange(xL+Dx/2.,xR,Dx)
    ngrid = len(grid)
    rho = scipy.ones_like(grid)/(xR-xL)

    # construction of sampler 
    param = Sampler.getDefaultParameters()
    sampler = Sampler(param)

    N = 500
    
    y = sampler.lift(rho,grid,N)
    dW = sampler.getBrownianIncrement(N)
    
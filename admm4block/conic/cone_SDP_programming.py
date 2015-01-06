import numpy as np

class ErrorDim(Exception):
    def __init__(self, dim1, dim2):
        self._dim1 = dim1
        self._dim2 = dim2
    def __str__(self):
        errstr = "Dimension should be "+repr(self._dim1)+" but is "+repr(self._dim2)
        return errstr

""" Defines a conic problem with:
     minimize Copt*x
     s.t. Aeq*x == beq
          Ain*x <= bin
          x is in K
          x is in Kp """
class ConicProgrammingProblem(object):

    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None, K=None, Kp=None):
        self._n = Copt.shape[0]
        self._neq = Aeq.shape[0]
        if Ain is not None:
            self._nin = Ain.shape[0]
        else:
            self._nin = 0
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        #K, Kp refer to PI K and PI Kp (projections), we use moreau's decomposition theorem
        self._Kp = Kp
        self._K = K
        
    def _Kpdual(self,z):
        return z + self._Kp(-z)
    def _Kdual(self,s,n):
        return s + self._K(-s,n)
        
    def __ADMM_noIn_step(self, X,s,z,y,AeqInv,sigma, tau):
        s = self._Kdual(self._Copt - z - np.dot(self._Aeq.T,y) - X/sigma,np.sqrt(self._n))
        y = np.dot(AeqInv,np.dot(self._Aeq,(self._Copt - s - z))) 
        z = self._Kpdual(self._Copt - s - np.dot(self._Aeq.T,y) - X/sigma)
        y = np.dot(AeqInv,np.dot(self._Aeq,(self._Copt - s - z))) 
        X += tau*sigma*(s + z + np.dot(self._Aeq.T,y) - self._Copt)
        return [X, s, z, y]

    def __ADMM_noIn(self, X0, s0, z0, AeqInv, sigma, tau,tol,nsteps):

        def CheckConditions(x,s,z,y):
            resP=np.sqrt(sum((np.dot(self._Aeq,x)-self._beq)**2))/(1+np.sqrt(sum(self._beq**2)))
            resD=np.sqrt(sum((s+z+np.dot(self._Aeq.T,y)-self._Copt)**2))/(1+np.sqrt(sum(self._Copt**2)))
            resK=np.sqrt(sum(x-self._K(x,np.sqrt(self._n))**2))/(1+np.sqrt(sum(x**2)))
            resKp=np.sqrt(sum((x-self._Kp(x))**2))/(1+np.sqrt(sum(x**2)))
            resKst=np.sqrt(sum(s-self._Kdual(s,np.sqrt(self._n))**2))/(1+np.sqrt(sum(s**2)))
            resKpst=np.sqrt(sum((z-self._Kpdual(z))**2))/(1+np.sqrt(sum(z**2)))
            resC1 = abs(np.dot(x,s))/(1+np.sqrt(sum(x**2))+np.sqrt(sum(s**2)))
            resC2 = abs(np.dot(x,z))/(1+np.sqrt(sum(x**2))+np.sqrt(sum(z**2)))
            return max(resP,resD,resK,resKp,resKst,resKpst,resC1,resC2)

        #tau should be less than (1+sqrt(5))/2 for convergence 
        y = np.dot(AeqInv, np.dot(self._Aeq,self._Copt - s0 - z0))
        res = CheckConditions(X0,s0,z0,y)
        k = 0
        while res > tol and k < nsteps:
            [X0, s0, z0, y] = self.__ADMM_noIn_step(X0,s0,z0,y,AeqInv,sigma,tau)
            res = CheckConditions(X0,s0,z0,y)
            k += 1
        return [X0, s0, z0, y, res]

    def InitConditions(self,x0=None,s0=None,z0=None): #Is this necessary?
        if x0 is None:
            x0 = np.dot(np.linalg.pinv(self._Aeq),self._beq)
        if s0 is None:
            s0 = x0
        if z0 is None:
            z0 = x0
        s0 = self._K(s0,np.sqrt(self._n))
        z0 = self._Kp(z0)
        return [x0,s0,z0]

    def Solve(self,sigma, tau, tol, nsteps,X0 = None, s0 = None, z0 = None, AeqInv = None):
        [X0,s0,z0] = self.InitConditions(X0,s0,z0)
        if AeqInv is None:
            AeqInv = np.linalg.inv(np.dot(self._Aeq,self._Aeq.T))
        if self._nin == 0:
            [x,s,z,y,res] = self.__ADMM_noIn(X0,s0,z0,AeqInv,sigma,tau,tol,nsteps)
            return [x,s,z,y,res,"no inequalities"]
        else:
            return "Not yet done"

    def Print(self, x, y):
        """Write the results of the primal and dual coefficients to an ouptput file. """
        f = open('output', 'w')
        f.write('The primal coefficients x are:')
        s_x = ''.join(str(e) for e in x) # assume x and y are lists
        f.write(s_x)
        f.write('The dual coefficients y are:')
        s_y = ''.join(str(e) for e in y)
        f.write(s_y)
        f.close()
            
"""Defines an SDP problem with:
     minimize Tr(Copt*X)
     s.t. tr(Aeq[k]*X) == beq[k] k = 0,...,M-1
          tr(Ain[k]*X) <= bin[k] k = 0,...,M'-1
          X is SDP
          Elements of X nonnegative"""
class DNNSDP(object):

    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None):
        self._n = Copt.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        
#number of columns of Aeq has to be n

    def toConic(self):
        if self._Aeq is not None:
            Aeq = self._Aeq[0].reshape(-1)
            for Mat in self._Aeq[1:]:
                Aeq = np.vstack((Aeq, Mat.reshape(-1)))
        
            beq = np.array(self._beq).transpose()
        else:
            Aeq = None
            beq = None
        if self._Ain is not None:
            Ain = self._Ain[0].reshape(-1)
            for Mat in self._Ain[1:]:
                Ain = np.vstack((Ain, Mat.reshape(-1)))

            bin = np.array(self._bin).transpose()
        else:
            Ain = None
            bin = None
        def Kp(X):
            X[X<0.] = 0
            return X

        def K(X,n):
            matX = X.reshape(n,n)
            B = np.linalg.eigh(matX)
            B[0][B[0]<0.] = 0.
            B[0][abs(B[0]<1e-9)] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T

        ConicP = ConicProgrammingProblem(self._Copt.reshape(-1),Aeq,beq,Ain,bin, K, Kp)
        return ConicP
        
    def Solve(self, sigma, tau, tol, nsteps,X0 = None, s0 = None, z0 = None, AeqInv = None):
        if X0 is not None:
            X0 = X0.reshape(-1)
        if s0 is not None:
            s0 = s0.reshape(-1)
        if z0 is not None:
            z0 = z0.reshape(-1)
        myCon = self.toConic()
        [X,s,z,y,res,mark] = myCon.Solve(sigma,tau,tol,nsteps,X0,s0,z0,AeqInv)
        return [X.reshape(self._n,self._n),s.reshape(self._n,self._n),z.reshape(self._n,self._n),y,res,mark]
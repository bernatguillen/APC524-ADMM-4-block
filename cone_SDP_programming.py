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
        #K, Kp refer to PI K* and PI Kp* (projections to the duals)
        self._Kp = Kp
        self._K = K

    def __ADMM_noIn_step(self, X,s,z,y,AeqInv,sigma, tau):
        s = self._K(self._Copt - z - np.dot(self._Aeq.T,y) - X/sigma,np.sqrt(self._n))
        y = np.dot(AeqInv,np.dot(self._Aeq,(self._Copt - s - z))) #does z change in K? Careful
        z = self._Kp(self._Copt - s - np.dot(self._Aeq.T,y) - X/sigma)
        y = np.dot(AeqInv,np.dot(self._Aeq,(self._Copt - s - z))) #does s change in Kp? Careful
        X += tau*sigma*(s + z + np.dot(self._Aeq.T,y) - self._Copt)
        return [X, s, z, y]

    def __ADMM_noIn(self, X0, s0, z0, AeqInv, sigma, tau,tol,nsteps):

        def CheckConditions(x,s,z,y):
            res1=abs(np.sqrt(sum((np.dot(self._Aeq,x)-self._beq)**2)))
            res2=abs(sum((s+z+np.dot(self._Aeq.T,y)-self._Copt)**2))
            res3=abs(sum(x*s))
            res4=abs(sum(x*z))
            return max(res1,res2,res3,res4)

        #tau should be less than (1+sqrt(5))/2 for convergence 
        y = np.dot(AeqInv, np.dot(self._Aeq,self._Copt - s0 - z0))
        res = CheckConditions(X0,s0,z0,y)
        k = 0
        while res > tol and k < nsteps:
            [X0, s0, z0, y] = self.__ADMM_noIn_step(X0,s0,z0,y,AeqInv,sigma,tau)
            res = CheckConditions(X0,s0,z0,y)
            k += 1
        if res <= tol:
            status = 0
        else:
            status = 1
        return [X0, s0, z0, y, status]

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
            [x,s,z,y,status] = self.__ADMM_noIn(X0,s0,z0,AeqInv,sigma,tau,tol,nsteps)
            return [x,s,z,y,status,"no inequalities"]
        else:
            return "Not yet done"
            
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
            matX = X.T.reshape(n,n)
            B = np.linalg.eigh(matX)
            B[0][B[0]<0.] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T

        ConicP = ConicProgrammingProblem(self._Copt.reshape(-1),Aeq,beq,Ain,bin, K, Kp)
        return ConicP
        """How do you want me to add Kp, K? As functions? Maybe add directly delta K* and delta K*?"""

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
class ConicProgramming(object):

    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None, K=None, Kp=None):
        self._n = Copt.shape[1]
        self._neq = Aeq.shape[0]
        self._nin = Ain.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        #K, Kp refer to PI K* and PI Kp* (projections to the duals)
        self._Kp = Kp
        self._K = K


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
        Aeq = self._Aeq[1].reshape(-1)
        for Mat in self._Aeq[2:]:
            Aeq = np.vstack((Aeq, Mat.reshape(self._n**2)))
        
        beq = np.array(self._beq).transpose()
        Ain = self._Ain[1].reshape(self._n**2)
        for Mat in self._Ain[2:]:
            Ain = np.vstack((Ain, Mat.reshape(self._n**2)))

        bin = np.array(self._bin).transpose()

        def Kp(X):
            X[X<0.] = 0

        def K(X,n):
            matX = X.T.reshape(n,n)
            B = np.linalg.eig(matX)
            B[0][B[0]<0.] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T

        ConicP = ConicProgramming(self._Copt.reshape(self._n**2),Aeq,beq,Ain,bin, K, Kp)
        """How do you want me to add Kp, K? As functions? Maybe add directly delta K* and delta K*?"""

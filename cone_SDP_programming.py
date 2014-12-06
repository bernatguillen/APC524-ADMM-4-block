import numpy as np

class ErrorDim(Exception):
    def __init__(self, dim1, dim2):
        self._dim1 = dim1
        self._dim2 = dim2
    def __str__(self):
        errstr = "Dimension should be "+repr(self._dim1)+" but is "+repr(self._dim2)
        return errstr


class ConicProgramming(object):
"""Defines a conic problem with:
     minimize Copt*x
     s.t. Aeq*x == beq
          Ain*x <= bin
          x is in K
          x is in Kp"""
    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None):
        self._n = Copt.shape[1]
        self._neq = Aeq.shape[0]
        self._nin = Ain.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        """How should I receive K and Kp here? Maybe receive functions delta K and delta Kp and create the functions delta K* and delta Kp*?"""


class DNNSDP(object):
"""Defines an SDP problem with:
     minimize Tr(Copt*X)
     s.t. tr(Aeq[k]*X) == beq[k] k = 0,...,M-1
          tr(Ain[k]*X) <= bin[k] k = 0,...,M'-1
          X is SDP
          Elements of X nonnegative"""
    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None):
        self._n = Copt.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        
#number of columns of Aeq has to be n
"""For now, no error checking on dimensions or number of constraints"""
    def toConic(self):
        Aeq = self._Aeq[1].reshape(self._n**2)
        for Mat in self._Aeq[2:]:
            Aeq = np.vstack((Aeq, Mat.reshape(self._n**2)))
        
        beq = np.matrix(self._beq).transpose()
        Ain = self._Ain[1].reshape(self._n**2)
        for Mat in self._Ain[2:]:
            Ain = np.vstack((Ain, Mat.reshape(self._n**2)))

        bin = np.matrix(self._bin).transpose()

        ConicP = ConicProgramming(self._Copt.reshape(self._n**2),Aeq,beq,Ain,bin)
        """How do you want me to add Kp, K? As functions? Maybe add directly delta K* and delta K*?"""

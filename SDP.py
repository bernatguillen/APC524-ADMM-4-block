import numpy as np

class ErrorDim(Exception):
    def __init__(self, dim1, dim2):
        self._dim1 = dim1
        self._dim2 = dim2
    def __str__(self):
        errstr = "Dimension should be "+repr(self._dim1)+" but is "+repr(self._dim2)
        return errstr

class SDPProblem(object):
    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None):
        self._n = Copt.shape[::-1]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin

        if self._Aeq.shape != self._n[::-1]:
            try:
                raise ErrorDim(self._n[::-1],self._Aeq.shape)
            except ErrorDim:
                print 'Error in equality matrix constraint'

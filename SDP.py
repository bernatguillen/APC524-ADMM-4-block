import numpy as np

class ErrorDim(Exception):
    def __init__(self, dim1, dim2):
        self._dim1 = dim1
        self._dim2 = dim2
    def __str__(self):
        errstr = "Dimension should be "+repr(self._dim1)+" but is "+repr(self._dim2)
        return errstr

class SDPProblem(object):
"""Defines an SDP problem with:
     minimize Tr(Copt*X)
     s.t. tr(Aeq[k]*X) == beq[k] k = 0,...,M-1
          tr(Ain[k]*X) <= bin[k] k = 0,...,M'-1
          X is SDP
          Elements of X may or may not be nonnegative"""
    def __init__(self, Copt=None, Aeq=None, beq=None, Ain=None, bin=None, doubly = False):
        self._n = Copt.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        self._doubly = doubly
#number of columns of Aeq has to be n
        if self._Aeq.shape[1] != self._n:
            try:
                raise ErrorDim(self._n,self._Aeq.shape)
            except ErrorDim:
                print 'Error in equality matrix constraint'
        if self._beq.shape[0] != self._Aeq.shape[0]:
            try:
                raise ErrorDim(self._beq.shape[0], self._Aeq.shape[0])
            except ErrorDim:
                print 'Error in equality matrix constraint'
        if self._Ain.shape[1] != self._n:
            try:
                raise ErrorDim(self._n, self._Ain.shape[0])
            except ErrorDim:
                print 'Error in inequality matrix constraint'
        if self._bin.shape[0] != self._Ain.shape[0]:
            try:
                raise ErrorDim(self._bin.shape[0], self._Ain.shape[0])
            except ErrorDim:
                print 'Error in inequality matrix constraint'

    def toConic(self):
        

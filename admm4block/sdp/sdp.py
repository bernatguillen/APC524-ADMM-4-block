# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 22:51:29 2015

@author: Bernat
"""
import numpy as np
from admm4block.conic.conic import ConicProgrammingProblem
from scipy import linalg

class ErrorEqs(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ErrorDim(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SDP(object):

    def __init__(self, Copt, Aeq, beq, Ain=None, bin=None):
        def Kp(X):
            return X

        def K(X,n):
            matX = X.reshape(n,n)
            B = linalg.eigh(matX)
            B[0][B[0]<0.] = 0.
            #B[0][abs(B[0]<1e-9)] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T
        
        if Copt.shape[0] != Copt.shape[1]:
            try:
                raise ErrorDim(1)
            except ErrorDim:
                print 'Copt should be a square matrix'
                return                
        if len(set([mat.shape for mat in Aeq])) != 1 or Aeq[0].shape != Copt.shape:
            try:
                raise ErrorDim(2)
            except ErrorDim:
                print 'There is a problem with the dimensions of the equality constraints'
                return
        if len(Aeq)!=len(beq):
            try:
                raise ErrorEqs(len(beq))
            except ErrorEqs:
                print 'len(Aeq) should be equal to len(beq)'
                return
        if Ain is not None and (len(set([mat.shape for mat in Ain])) != 1 or Ain[0].shape != Copt.shape):
            try:
                raise ErrorDim(2)
            except ErrorDim:
                print 'There is a problem with the dimensions of the inequality constraints'
                return
        if Ain is not None and len(Ain)!=len(bin):
            try:
                raise ErrorEqs(len(beq))
            except ErrorEqs:
                print 'len(Ain) should be equal to len(bin)'
                return

        self._n = Copt.shape[0]
        self._Copt = Copt
        self._Aeq = Aeq
        self._beq = beq
        self._Ain = Ain
        self._bin = bin
        self._K = K
        self._Kp = Kp
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

        ConicP = ConicProgrammingProblem(self._Copt.reshape(-1),Aeq,beq, self._K, self._Kp,Ain,bin)
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
        
class DNNSDP(SDP):
    def __init__(self,Copt, Aeq, beq, Ain=None, bin=None):
        SDP.__init__(self,Copt, Aeq, beq, Ain, bin)
        def Kp(X):
            X[X<0.]=0.
            return X
        self._Kp = Kp
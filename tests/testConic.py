# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 11:09:26 2015

@author: Bernat
"""

import unittest

import numpy as np
import admm4block

class TestConic(unittest.TestCase):
    def testConeProjection(self):
        def K(X,n):
            B = np.linalg.eigh(X)
            B[0][B[0]<0.] = 0.
            #B[0][abs(B[0]<1e-9)] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C
        X0 = np.array([[-1.,0.],[0.,1.]])
        np.testing.assert_array_almost_equal(K(X0,3),np.array([[0.,0.],[0.,1.]]))
        
    def testSelfDualCone(self):
        def K(X):
            X[X<0.] = 0.
            return X
        def Kdual(X):
            return X + K(-X)
        X0 = np.array([-1., 2., 3.])
        np.testing.assert_array_almost_equal(K(X0),Kdual(X0))
        
    def testStep(self):
        def Kp(X):
            X[X<0.] = 0.
            return X
        def K(X,n):
            matX = X.reshape(n,n)
            B = np.linalg.eigh(matX)
            B[0][B[0]<0.] = 0.
            #B[0][abs(B[0]<1e-9)] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T
        
        X0 = np.array([1.,2.,2.,2.])
        s0 = K(X0,2)
        X = admm4block.ConicProgrammingProblem(Copt = np.array([1.,1.,1.,1.]), Aeq = np.array([[1.,0.,0.,0.]]), beq = np.array([1.]), K= K, Kp=Kp)
        [Xf,_,_,_,_,_]=X.Solve(sigma=1,tau=1,tol=1,nsteps = 1, X0=X0, s0=s0, z0 = X0)
        np.testing.assert_array_almost_equal(Xf , np.array([1.,0.07873759,0.07873759,1.54084855]))
    def testSolvelown(self):
        def Kp(X):
            X[X<0.] = 0.
            return X
        def K(X,n):
            matX = X.reshape(n,n)
            B = np.linalg.eigh(matX)
            B[0][B[0]<0.] = 0.
            #B[0][abs(B[0]<1e-9)] = 0.
            C = np.dot(B[1], (B[0]*B[1]).T)
            return C.reshape(-1).T
        
        X0 = np.array([1.,2.,2.,2.])
        s0 = K(X0,2)
        X = admm4block.ConicProgrammingProblem(Copt = np.array([1.,1.,1.,1.]), Aeq = np.array([[1.,0.,0.,0.]]), beq = np.array([1.]), K= K, Kp=Kp)
        [Xf,_,_,_,_,_]=X.Solve(sigma=1,tau=1,tol=1e-6,nsteps = 100, X0=X0, s0=s0, z0 = X0)
        np.testing.assert_array_almost_equal(Xf , np.array([1.,0.00,0.00,0.00]))

if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 10 11:09:26 2015

@author: Bernat
"""

import unittest

import numpy as np
from scipy import linalg
import admm4block

class TestSDP(unittest.TestCase):
        
    def testStepDNN(self):
        X0 = np.array([[1.,2.],[2.,2.]])
        X = admm4block.DNNSDP(Copt = np.array([[1.,1.],[1.,1.]]), Aeq = [np.array([[1.,0.],[0.,0.]])], beq = [1.])
        [Xf,_,_,_,_,_]=X.Solve(sigma=1,tau=1,tol=0.1,nsteps = 1, X0=X0,AeqInv=1.)
        np.testing.assert_array_almost_equal(Xf , np.array([[1.,0.07873759],[0.07873759,1.54084855]]))
    def testSolvelownDNN(self):
        X0 = np.array([[1.,2.],[2.,2.]])
        X = admm4block.DNNSDP(Copt = np.array([[1.,1.],[1.,1.]]), Aeq = [np.array([[1.,0.],[0.,0.]])], beq = [1.])
        [Xf,_,_,_,_,_]=X.Solve(sigma=1,tau=1,tol=0.1,nsteps = 100, X0=X0,AeqInv=1.)
        np.testing.assert_array_almost_equal(Xf , np.array([[1.,0.],[0.,0.]]))
    def testDNNpositivity(self):
        X0 = np.array([[1.,2.],[2.,2.]])
        X = admm4block.DNNSDP(Copt = np.array([[1.,1.],[1.,1.]]), Aeq = [np.array([[1.,0.],[0.,0.]])], beq = [1.])
        [Xf,_,_,_,_,_]=X.Solve(sigma=1,tau=1,tol=0.1,nsteps = 1, X0=X0)
        np.testing.assert_array_less(0,Xf)
    def testSDP(self):
        Aeq = [np.array([[1.,0.,1.],[0.,3.,7.],[1.,7.,5.]]),np.array([[0.,2.,8.],[2.,6.,0.],[8.,0.,4.]])]
        beq = [11.,19.]
        Copt = np.array([[1.,2.,3.],[2.,9.,0.],[3.,0.,7.]])
        mySDP = admm4block.SDP(Copt,Aeq,beq)
        [Xf,_,_,_,_,_] = mySDP.Solve(2,2,0.1,100)
        np.testing.assert_array_almost_equal(mySDP._K(Xf.reshape(-1),3).reshape(3,3),Xf)
        
if __name__ == '__main__':
    unittest.main()
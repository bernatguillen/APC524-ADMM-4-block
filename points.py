# -*- coding: utf-8 -*-
"""
Created on Tue Jan 06 02:17:46 2015

@author: Bernat
"""
#%%

import numpy as np
import admm4block

def f(x1,x2,s):
    return np.exp(-np.sqrt(sum((x1-x2)**2))/s)

def W(points,s):
    A = np.ndarray((len(points),len(points)))
    for i in range(len(points)):
        for j in range(len(points)):
            A[i][j] = f(points[i],points[j],s)
    return A

points = [np.array([0.,0.]),np.array([6.,6.]),np.array([0.,1.]),np.array([5.,5.]),np.array([4.,5.]),np.array([5.,4.])]
Copt= -W(points,4)
Aeq = [np.zeros((6,6)) for i in range(6)] + [np.eye(6)]
for i in range(6):
    Aeq[i][:,i] = 1.
beq = [1.,1.,1.,1.,1.,1.,2.]

mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(10,0.6,0.1,10000)
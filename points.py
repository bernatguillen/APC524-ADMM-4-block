# -*- coding: utf-8 -*-
"""
Created on Tue Jan 06 02:17:46 2015

@author: Bernat
"""
#%%

import numpy as np
import admm4block
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def f(x1,x2,s):
    return np.exp(-np.sqrt(sum((x1-x2)**2))/s)

def W(points,s):
    A = np.ndarray((len(points),len(points)))
    for i in range(len(points)):
        for j in range(len(points)):
            A[i][j] = f(points[i],points[j],s)
    return A
k = 8
n = 24
means = [4*np.array([np.cos(i*2*np.pi/k),np.sin(i*2*np.pi/k)]) for i in range(k)]
covs =  [np.eye(2)/2. for i in range(k)]
for i in range(k):
    covs[i][0][1] = covs[i][1][0] = np.random.rand(1)[0]
points = [np.random.multivariate_normal(means[i],covs[i],n/k) for i in range(k)]
points = [vec for sublist in points for vec in sublist]


Copt= -W(points,4)
Aeq = [np.zeros((n,n)) for i in range(n)] + [np.eye(n)]
for i in range(n):
    Aeq[i][:,i] = 1.
beq = [1.]*n + [k]

mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(1,1,0.1,1000)

points = np.array(points)
plt.scatter(points[:,0],points[:,1])

colors = iter(cm.rainbow(np.linspace(0, 1, k)))
s = set()
for i in range(n):
    if i not in s:
        idx = np.where(X[:,i]>=1./n)[0].tolist()
        s = s.union(set(idx))
        plt.scatter(points[idx,0],points[idx,1],color = next(colors))
        

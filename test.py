#!/usr/bin/env python
import admm4block
from make_1Dsig import make_1Dsig
import numpy as np
import math

L = 5
time = 3
noise = 0.0
sigma = 1
nsteps = 1000
tau = 0.1
tol = 1e-3
(movie,shifts) = make_1Dsig(L,time,noise)
Copt = np.reshape(movie,(L*time,1),'F')
Copt = np.multiply(Copt,Copt.T)
neq = time*1/2*L*(L+1) + (1/2*time*(time-1))*(L+L*(L-1))
Aeq = np.random.randn(neq,L*L*time*time)
Aeq = [Aeq[i,:].reshape((L*time,L*time)) for i in range(neq)]
beq = np.random.randn(1,neq).tolist()[0]
#print np.dot(np.linalg.pinv(Aeq),beq.T).shape
mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(sigma, tau, tol, nsteps)

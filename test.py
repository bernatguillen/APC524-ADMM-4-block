#!/usr/bin/env python
import admm4block
from make_1Dsig import make_1Dsig
import numpy as np
import math

L = 3
time = 3
noise = 0.0
sigma = 1
nsteps = 1000
tau = 0.1
tol = 1e-3
def index(i,k):
    return i + k*L
(movie,shifts) = make_1Dsig(L,time,noise)
Copt = np.reshape(movie,(L*time,1),'F')
Copt = np.multiply(Copt,Copt.T)
Aeq = []
beq = []
for t in range(0,time):
    for i in range(0,L):
        for j in range(i,L):
            A = np.zeros((L*time,L*time))
            A[index(i,t),index(j,t)] = 1
            Aeq.append(A)
            if (i == j):
                beq.append(1)
            else:
                beq.append(0)
for k in range(0,time-1):
    for l in range(k+1,time):
        for i in range(0,L):
            A = np.zeros((L*time,L*time))
            A[index(i,k),index(0,l):index(L,l)] = 1;
            Aeq.append(A)
            beq.append(1)
        for i in range(0,L):
            for j in range(0,L-1):
                A = np.zeros((L*time,L*time))
                A[index((i+j)%L,k),index(j,l)] = 1
                A[index((i+j+1)%L,k),index(j+1,l)] = -1
                Aeq.append(A)
                beq.append(0)

mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(sigma, tau, tol, nsteps)

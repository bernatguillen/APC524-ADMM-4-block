#!/usr/bin/env python
import admm4block
import numpy as np
import math
import time as timing
from mlab.releases import latest_release

def make_1Dsig(time,L,sigma):
    I = np.random.randn(L,1)
    I = math.sqrt(L)*I/np.linalg.norm(I)
    movie = np.zeros((L,time))
    shifts = np.random.randint(L, size=time)
    for t in range(0,time):
        shift = shifts[t]
        movie[:,t:t+1] = np.roll(I,-shift) + sigma*np.random.randn(L,1)
    return (movie,shifts)

# Initial parameters for the multireference alignment SDP problem
L = 5
time = 5
noise = 0.0
def index(i,k):
    return i + k*L
(movie,shifts) = make_1Dsig(L,time,noise)

# Solve using our ADMM package
# First build Copy, Aeq, and beq matrices
Copt = np.reshape(movie,(L*time,1),'F')
Copt = np.multiply(Copt,Copt.T)
start_admm = timing.clock()
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
	for i in range(1):
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
# Conditions for the solver
sigma = 1
tau = 1
tol = 1e-3
nsteps = 1000
mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(sigma, tau, tol, nsteps)
end_admm = timing.clock()
print "Time used by ADMM: " + str(end_admm - start_admm) + " s"

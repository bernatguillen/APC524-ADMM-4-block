#!/usr/bin/env python
import admm4block
import numpy as np
import math

def make_1Dsig(time,L,sigma):
    I = np.random.randn(L,1)
    I = math.sqrt(L)*I/np.linalg.norm(I)
    movie = np.zeros((L,time))
    shifts = np.random.randint(L, size=time)
    for t in range(0,time):
        shift = shifts[t]
        movie[:,t:t+1] = np.roll(I,-shift) + sigma*np.random.randn(L,1)
    return (movie,shifts)
    
L = 2
time = 2
noise = 0.0
sigma = 1
nsteps = 1000
tau = 1
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
	for i in range(1):
        #for i in range(0,L):
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
#Aeq1 = Aeq[0].reshape(-1)
#for Mat in Aeq[1:]:
#    Aeq1 = np.vstack((Aeq1, Mat.reshape(-1)))

#beq = np.array(beq).transpose()
#print np.dot(np.linalg.pinv(Aeq1),beq).reshape(L*time,L*time)
[X,s,z,y,res,_]=mySDP.Solve(sigma, tau, tol, nsteps)
print X

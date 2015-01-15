#!/usr/bin/env python
import admm4block
import numpy as np
import math
import time as timing
import warnings
from mlabwrap import mlab

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
L = 2
time = 2
noise = 0.5
def index(i,k):
    return i + k*L
(movie,shifts) = make_1Dsig(L,time,noise)

########################### Solve using our ADMM package ###########################
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
            A[index(i,k),index(0,l):index(L,l)] = 1
            Aeq.append(A)
            beq.append(1)
        for i in range(0,L):
            for j in range(0,L-1):
                A = np.zeros((L*time,L*time))
                A[index((i+j)%L,k),index(j,l)] = 1
                A[index((i+j+1)%L,k),index(j+1,l)] = -1
                Aeq.append(A)
                beq.append(0)
#For debugging purposes
#for k in range(0,time-1):
#	for l in range(k+1,time):
#		for i in range(0,L):
#			for j in range(0,L):
#				A = np.zeros((L*time,L*time))
#				A[index(i,k),index(j,l)] = 1
#				A[index(j,l),index(i,k)] = -1
#				Aeq.append(A)
#				beq.append(0)

# Conditions for the solver
sigma = 1
tau = 1
# See http://cvxr.com/cvx/doc/solver.html for more information on CVX precision
machine_precision = 2.2e-16
tol_low = machine_precision**0.25       # Call alignment.m with 1 as second argument
tol_default = machine_precision**0.375  # Call alignment.m with 2 as second argument
tol_high = machine_precision**0.5       # Call alignment.m with 3 as second argument
nsteps = 1000
start_solve = timing.clock()
mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
[X,s,z,y,res,_]=mySDP.Solve(sigma, tau, tol_low, nsteps)
end_admm = timing.clock()
# Time to solve the SDP using ADMM (already given Copt, Aeq, and beq)
print "Solve time for ADMM: " + str(end_admm - start_solve) + " s"
# Total time (includes solving using ADMM and building the Copt, Aeq, and beq matrices)
print "Total time for ADMM: " + str(end_admm - start_admm) + " s"
print ""
################################# Solve using CVX #################################

# Suppress ComplexWarning so that output is not messed up
warnings.filterwarnings("ignore")
Xcvx,cvxtime = mlab.alignment(movie,1,nout=2)
cvxtime = np.round(cvxtime, 6)
# Only solving
print "Solve time for CVX : " + str(cvxtime[0,4]) + " s"
# Total time (includes model building, presolving, and solving)
print "Total time for CVX : " + str(sum(cvxtime[0,2:5])) + " s"
#print X
#print Xcvx

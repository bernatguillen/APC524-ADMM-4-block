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


import numpy as np

def pbc1(x,L):
    if x > L/2: x -= L
    elif x<-L/2: x += L

    return x

def pbc2(x,L):
    x = x - np.rint(x/L)*L

    return x


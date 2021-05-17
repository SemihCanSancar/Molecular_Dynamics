import numpy as np

## putting the atoms array back into the simulation box
def pbc(x,L):
    x = x - np.rint(x/L)*L

    return x

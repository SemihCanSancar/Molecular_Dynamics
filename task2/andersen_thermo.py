import numpy as np

def anderson_thermo(vel,nu,temp):
    sigma = np.sqrt(temp)

    for i in range(len(vel)):
        if np.random.ranf() < nu:
            
            vel[i] = np.random.normal(0,sigma,3)

    return vel
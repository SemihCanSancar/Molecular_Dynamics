import numpy as np

def average(energy_array,start):
    energy_array = energy_array[start:]
    return np.sum(energy_array)/len(energy_array)

def variance(avg,energy_array,start):
    energy_array = energy_array[start:]
    l = len(energy_array)
    erg = 0
    for i in range(l):
        erg += (energy_array[i]-avg)**2
    
    return np.sqrt(1/(l-1)*erg)/np.sqrt(l)


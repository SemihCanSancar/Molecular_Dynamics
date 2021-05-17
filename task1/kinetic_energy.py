import numpy as np

def kinetic_en(vel_arr):

    kin_en = 0.5*np.sum(vel_arr**2)

    return kin_en
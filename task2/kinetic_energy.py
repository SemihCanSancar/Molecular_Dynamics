import numpy as np

def kinetic_en(vel_arr):
    kin_en = 0

    for i in range(len(vel_arr)):
        kin_en += 0.5*(vel_arr[i][0]**2 + vel_arr[i][1]**2 + vel_arr[i][2]**2)

    return kin_en
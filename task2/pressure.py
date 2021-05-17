import numpy as np

def pressure(density,box_len,kT,atoms_arr,forces_arr):

    p1 = kT*density
    p2 = 0

    for i in range(len(atoms_arr)):
        for j in range(len(atoms_arr)):
            if i != j:
                p2 += atoms_arr[i][0]*forces_arr[j][0]+atoms_arr[i][1]*forces_arr[j][1]+atoms_arr[i][2]*forces_arr[j][2]

    return p1+p2
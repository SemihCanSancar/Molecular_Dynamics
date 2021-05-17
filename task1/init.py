import numpy as np
import pbc 

## initialize simple cube geometry

def init_position(n_atoms,density):
    atoms_arr = np.zeros((n_atoms,3))
    num_atoms_dim = round((n_atoms)**(1/3))
    box_len = num_atoms_dim/(density**(1/3))
    dist = box_len/num_atoms_dim
    x_vec = np.transpose([1,0,0])
    y_vec = np.transpose([0,1,0])
    z_vec = np.transpose([0,0,1])

    i = 0
    for x in range(num_atoms_dim):
        for y in range(num_atoms_dim):
            for z in range(num_atoms_dim):
                atoms_arr[i] = x*x_vec*dist + y*y_vec*dist + z*z_vec*dist  
                atoms_arr[i] = atoms_arr[i]
                i += 1


    ## putting the cube into the simulation box
    atoms_arr = pbc.pbc(atoms_arr,box_len)
    
    return atoms_arr,box_len
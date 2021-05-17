import numpy as np

def volume(r,dr):
    return (4*np.pi/3)*((r+dr/2)**3-(r-dr/2)**3)
    

def number_particles(atoms_arr,r,dr,box_len):
    
    ## conditions from the middle point of the box
    state1 = (box_len/2) - atoms_arr[:,0:1] > r
    state2 = (box_len/2) - atoms_arr[:,0:1] < (box_len - (dr+r)) 
    state3 = (box_len/2) - atoms_arr[:,0:1] > r
    state4 = (box_len/2) - atoms_arr[:,0:1] < (box_len - (dr+r))
    state5 = (box_len/2) - atoms_arr[:,0:1] > r
    state6 = (box_len/2) - atoms_arr[:,0:1] < (box_len - (dr+r))

    # return indices in array that fulfil the conditions with numpy where
    indices = np.where(state1 * state2 * state3 * state4 * state5 * state6)
    
    num_atoms = len(indices[0])
    #print(num_atoms)
    return num_atoms


def g(r,dr,density,atoms_arr,box_len):
    erg = []
    for x in range(10):
        num_particles = number_particles(atoms_arr,dr,r,box_len)
        g = num_particles/(density*volume(r,dr))
        r += dr
        erg.append(g)
    
    return np.array(erg)
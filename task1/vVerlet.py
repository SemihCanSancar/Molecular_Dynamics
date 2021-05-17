import numpy as np
import lennard_jones as flj

def velocity_verlet(atoms_arr,vel_arr,forces_arr,dt,box_len):
    
    atoms_arr = atoms_arr + vel_arr*dt + 0.5*forces_arr*dt*dt
    vel_arr = vel_arr + forces_arr*dt*0.5
    forces_arr,pot = flj.force_LJ(atoms_arr, box_len)
    vel_arr = vel_arr + forces_arr*dt*0.5
    
    return atoms_arr,vel_arr,forces_arr, pot
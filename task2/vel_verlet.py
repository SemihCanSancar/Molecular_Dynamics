import numpy as np
import forces_LJ as flj

def velocity_verlet_andersen_thermo(atoms_arr,vel_arr,forces_arr,dt,box_len,nu,temp):
    
    def anderson_thermo(vel,nu,temp):
        sigma = np.sqrt(temp)

        for i in range(len(vel)):
            if np.random.ranf() < nu*dt:
            
                vel[i] = np.random.normal(0,sigma,3)

        return vel

    atoms_arr = atoms_arr + vel_arr*dt + 0.5*forces_arr*dt*dt
    vel_arr = vel_arr + forces_arr*dt*0.5
    forces_arr,pot = flj.force_LJ(atoms_arr, box_len)
    vel_arr = vel_arr + forces_arr*dt*0.5
    vel_arr = anderson_thermo(vel_arr,nu,temp)

    
    return atoms_arr,vel_arr,forces_arr, pot
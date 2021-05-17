import numpy as np
import forces_LJ as flj
import pbc
import scipy.constants as sc
import init_vel as iv
import init_pos as ip
import vel_verlet as vv

## melting system
def melt(atoms_arr,temp,box_len):
    vel_arr = iv.initialize(len(atoms_arr),temp)
    forces_arr,_ = flj.force_LJ(atoms_arr, box_len)
    for i in range(100):
        atoms_arr ,_,forces_arr,_ =  vv.velocity_verlet_andersen_thermo(atoms_arr,vel_arr,forces_arr,0.001,box_len,0.1,temp)
        atoms_arr = pbc.pbc2(atoms_arr,box_len)
    return atoms_arr

def write_pos(atoms_arr):
    trajec = open("init_iso_syst_i.xyz",'w')
    for t in range(1):
        trajec.write('%s \n \n' % (len(atoms_arr)))
        for i in range(len(atoms_arr)):
            trajec.write("Xe %s %s %s \n" % (atoms_arr[i][0],atoms_arr[i][1],atoms_arr[i][2]))

    trajec.close()


if __name__ == "__main__":
    
    density = [0.1,0.2,0.4,0.6,0.8]
    temp = 1.5

    ## melting the crystal structure
    high_temp = 1000


    ## initializing configuration of system 
    for dens in density:
        atoms_arr, box_len = ip.init_position(125, dens)

        atoms_arr = melt(atoms_arr,high_temp,box_len)

        np.save("atoms_arr_%.1f" %(dens),atoms_arr)
        np.save("box_len_%.1f" %(dens),box_len)
    
import vel_verlet as vv
import pbc as pbc
import forces_LJ as flj
import init_pos as ip
import init_vel as iv
import kinetic_energy as ke
import pressure as p
import radial_dist as rd
import average_variance as av

import numpy as np
import matplotlib.pyplot as plt

def simulation(timesteps, atoms_arr, vel_arr,box_len,nu,temp,density):
    pot_energy = []
    kin_energy = []
    pressure = []

    forces_arr,_ = flj.force_LJ(atoms_arr, box_len)

    for i in range(timesteps):

        # calculating position array, velocity array and potential energy
        atoms_arr, vel_arr,forces_arr, pot = vv.velocity_verlet_andersen_thermo(atoms_arr,vel_arr,forces_arr,0.00001,box_len,nu,temp)
        atoms_arr = pbc.pbc2(atoms_arr,box_len)

        pressure.append(p.pressure(density,box_len,temp,atoms_arr,forces_arr))

        pot_energy.append(pot)

        kin = ke.kinetic_en(vel_arr)
        kin_energy.append(kin)

    return atoms_arr,vel_arr,pot_energy,kin_energy, pressure

if __name__ == "__main__":
    
    ## parameters

    densities = [0.1,0.2,0.4,0.6,0.8]
    temp = 1.5

    ## empty energy lists to plot it at last and see density dependency

    pot_en = []
    kin_en = []
    total_en = []
    pressure = []
    radial = []
    ## loading the positions from previous step

    for dens in densities:
        ## reading the position arrays from the last point of melting the liquid crystal
        ## Length of the box loaded as well and velocities initialized
        atoms_arr = np.load("atoms_arr_%0.1f.npy" %(dens),mmap_mode="r+")
        box_len = np.load("box_len_%0.1f.npy" % (dens),mmap_mode="r+" )
        vel_arr = iv.initialize(len(atoms_arr),temp)

        print("Working on density: %0.1f ... \n \n" % (dens))

        atoms_arr,_,pot,kine,pres = simulation(1000,atoms_arr,vel_arr,box_len,0.1,temp,dens)
        
        ## energies
        pot = np.array(pot)
        pot_en.append(pot)
        kine = np.array(kine)
        kin_en.append(kine)
        total = pot + kine
        total_en.append(total)

        ##pressure
        pres = np.array(pres)
        pressure.append(pres)

        ## radial distribution
        radial.append(rd.g(1,0.25,dens,atoms_arr,box_len))

        
    ## averages and variance calculation
    len_en_arr = int(len(kin_en[0])/2)
    avg_kin_01 = av.average(kin_en[0],len_en_arr)
    avg_kin_02 = av.average(kin_en[1],len_en_arr)
    avg_kin_03 = av.average(kin_en[2],len_en_arr)
    avg_kin_04 = av.average(kin_en[3],len_en_arr)
    avg_kin_05 = av.average(kin_en[4],len_en_arr)

    var_kin_01 = av.variance(avg_kin_01,kin_en[0],len_en_arr)
    var_kin_02 = av.variance(avg_kin_02,kin_en[1],len_en_arr)
    var_kin_03 = av.variance(avg_kin_03,kin_en[2],len_en_arr)
    var_kin_04 = av.variance(avg_kin_04,kin_en[3],len_en_arr)
    var_kin_05 = av.variance(avg_kin_05,kin_en[4],len_en_arr)

    print("Kinetic energy with density 0.1: %.3f +/- %.3f" %(avg_kin_01,var_kin_01))
    print("Kinetic energy with density 0.2: %.3f +/- %.3f" %(avg_kin_02,var_kin_02))
    print("Kinetic energy with density 0.4: %.3f +/- %.3f" %(avg_kin_03,var_kin_03))
    print("Kinetic energy with density 0.6: %.3f +/- %.3f" %(avg_kin_04,var_kin_04))
    print("Kinetic energy with density 0.8: %.3f +/- %.3f \n \n" %(avg_kin_05,var_kin_05))
    
    
    avg_pot_01 = av.average(pot_en[0],len_en_arr)
    avg_pot_02 = av.average(pot_en[1],len_en_arr)
    avg_pot_03 = av.average(pot_en[2],len_en_arr)
    avg_pot_04 = av.average(pot_en[3],len_en_arr)
    avg_pot_05 = av.average(pot_en[4],len_en_arr)

    var_pot_01 = av.variance(avg_pot_01,pot_en[0],len_en_arr)
    var_pot_02 = av.variance(avg_pot_02,pot_en[1],len_en_arr)
    var_pot_03 = av.variance(avg_pot_03,pot_en[2],len_en_arr)
    var_pot_04 = av.variance(avg_pot_04,pot_en[3],len_en_arr)
    var_pot_05 = av.variance(avg_pot_05,pot_en[4],len_en_arr)


    print("Potential energy with density 0.1: %.3f +/- %.3f" %(avg_pot_01,var_pot_01))
    print("Potential energy with density 0.2: %.3f +/- %.3f" %(avg_pot_02,var_pot_02))
    print("Potential energy with density 0.4: %.3f +/- %.3f" %(avg_pot_03,var_pot_03))
    print("Potential energy with density 0.6: %.3f +/- %.3f" %(avg_pot_04,var_pot_04))
    print("Potential energy with density 0.8: %.3f +/- %.3f \n \n" %(avg_pot_05,var_pot_05))

    avg_tot_01 = av.average(total_en[0],len_en_arr)
    avg_tot_02 = av.average(total_en[1],len_en_arr)
    avg_tot_03 = av.average(total_en[2],len_en_arr)
    avg_tot_04 = av.average(total_en[3],len_en_arr)
    avg_tot_05 = av.average(total_en[4],len_en_arr)

    var_tot_01 = av.variance(avg_tot_01,total_en[0],len_en_arr)
    var_tot_02 = av.variance(avg_tot_02,total_en[1],len_en_arr)
    var_tot_03 = av.variance(avg_tot_03,total_en[2],len_en_arr)
    var_tot_04 = av.variance(avg_tot_04,total_en[3],len_en_arr)
    var_tot_05 = av.variance(avg_tot_05,total_en[4],len_en_arr)

    print("Total energy with density 0.1: %.3f +/- %.3f" %(avg_tot_01,var_tot_01))
    print("Total energy with density 0.2: %.3f +/- %.3f" %(avg_tot_02,var_tot_02))
    print("Total energy with density 0.4: %.3f +/- %.3f" %(avg_tot_03,var_tot_03))
    print("Total energy with density 0.6: %.3f +/- %.3f" %(avg_tot_04,var_tot_04))
    print("Total energy with density 0.8: %.3f +/- %.3f \n \n" %(avg_tot_05,var_tot_05))

    ## plot the energies
    
    plt.figure(1)
    plt.xlabel("timesteps")
    plt.ylabel("potential energy [kJ/mol]")
    plt.title("Potential energy plot")
    plt.plot(pot_en[0],color="red",label="density: 0.1 g/cm^3")
    plt.plot(pot_en[1],color="green",label="density: 0.2 g/cm^3")
    plt.plot(pot_en[2],color="blue",label="density: 0.4 g/cm^3")
    plt.plot(pot_en[3],color="black",label="density: 0.6 g/cm^3")
    plt.plot(pot_en[4],color="gray",label="density: 0.8 g/cm^3")
    plt.legend()
    plt.savefig("potential_energy_test.png")

    plt.figure(2)
    plt.xlabel("timesteps")
    plt.ylabel("kinetic energy [kJ/mol]")
    plt.title("Kinetic energy plot")
    plt.plot(kin_en[0],color ="red",label="density: 0.1 g/cm^3")
    plt.plot(kin_en[1],color ="green",label="density: 0.2 g/cm^3")
    plt.plot(kin_en[2],color ="blue",label="density: 0.4 g/cm^3")
    plt.plot(kin_en[3],color ="black",label="density: 0.6 g/cm^3")
    plt.plot(kin_en[4],color ="gray",label="density: 0.8 g/cm^3")
    plt.legend()
    plt.savefig("kinetic_energy_test.png")

    plt.figure(3)
    plt.xlabel("timesteps")
    plt.ylabel("total energy [kJ/mol]")
    plt.title("Total energy plot")
    plt.plot(total_en[0],color = "red",label="density: 0.1 g/cm^3")
    plt.plot(total_en[1],color = "green",label="density: 0.2 g/cm^3")
    plt.plot(total_en[2],color = "blue",label="density: 0.4 g/cm^3")
    plt.plot(total_en[3],color = "black",label="density: 0.6 g/cm^3")
    plt.plot(total_en[4],color = "gray",label="density: 0.8 g/cm^3")
    plt.legend()
    plt.savefig("total_energy_test.png")   

    ##plot pressure
    plt.figure(4)
    plt.xlabel("Timesteps")
    plt.ylabel("Pressure [Pascal]")
    plt.title("Pressure of the systems")
    plt.plot(pressure[0],color = "red",label="density: 0.1 g/cm^3")
    plt.plot(pressure[1],color = "green",label="density: 0.2 g/cm^3")
    plt.plot(pressure[2],color = "blue",label="density: 0.4 g/cm^3")
    plt.plot(pressure[3],color = "black",label="density: 0.6 g/cm^3")
    plt.plot(pressure[4],color = "gray",label="density: 0.8 g/cm^3")
    plt.legend()
    plt.savefig("pressure_test.png") 
    
    x = np.linspace(1,10,10)*0.25

    plt.figure(5)
    plt.plot(x,radial[0],color = "red",label="density: 0.1 g/cm^3")
    plt.plot(x,radial[1],color = "green",label="density: 0.2 g/cm^3")
    plt.plot(x,radial[2],color = "blue",label="density: 0.4 g/cm^3")
    plt.plot(x,radial[3],color = "black",label="density: 0.6 g/cm^3")
    plt.plot(x,radial[4],color = "gray",label="density: 0.8 g/cm^3")
    plt.xlabel(r"r [$\AA$]")
    plt.ylabel("g(r)")
    plt.title("Radial distribution function")
    plt.legend()
    plt.savefig("radial_dist_test_test.png")
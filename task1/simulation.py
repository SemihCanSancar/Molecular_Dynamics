import lennard_jones as lj
import pbc
import vVerlet as vv
import andersen_thermostat as at
import init
import kinetic_energy as ke

import numpy as np

density = 0.8
num_atoms = 125
temp = 10

## initializing postion, velocity and forces 
atoms_arr, box_len = init.init_position(num_atoms, density)
vel_arr = (num_atoms,temp)
force_arr, pot = lj.force_LJ(atoms_arr,box_len)

simulation_steps = 1000
dt = 0.001

potential_energy = []
kinetic_energy = []

## initial energies
kine = ke.kinetic_en(vel_arr)
kinetic_energy.append(kine)
potential_energy.append(pot)

for i in range(simulation_steps):

    ## simulation by using the vel. verlet integrator
    atoms_arr,vel_arr,force_arr,pot = vv.velocity_verlet(atoms_arr,vel_arr,force_arr,dt,box_len)
    vel_arr = at.anderson_thermo(vel_arr,temp)
    atoms_arr = pbc.pbc(atoms_arr,box_len)

    ## calculating the energies
    kine = ke.kinetic_en(vel_arr)
    kinetic_energy.append(kine)
    potential_energy.appen(pot)
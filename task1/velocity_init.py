import numpy as np

def initialize(n_particles,temp):
    sigma = np.sqrt(temp)

    vx = np.random.normal(0.0,sigma,n_particles)
    vy = np.random.normal(0.0,sigma,n_particles)
    vz = np.random.normal(0.0,sigma,n_particles)

    ## center of mass not drifting

    vx -= sum(vx)/n_particles
    vy -= sum(vx)/n_particles
    vz -= sum(vx)/n_particles

    ## rescaling the temp

    scale = np.sqrt(3*temp*n_particles/sum(vx**2+vy**2+vz**2))
    vx *= scale
    vy *= scale
    vz *= scale

    ## put velocities together as one matrix
    velocity = np.transpose((vx,vy,vz))
    
    return velocity
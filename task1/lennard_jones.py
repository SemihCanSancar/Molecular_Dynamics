import numpy as np
import pbc 

def force_LJ(r,box_len):
    N = len(r)
    force = np.zeros((N,3))
    cutoff = box_len/2
    pot = 0
    
    for i in range(N):
        for j in range(N):
            
            if i != j:
                dx = r[i][0] - r[j][0]
                dx = pbc.pbc1(dx,box_len)
                dy = r[i][1] - r[j][1]
                dy = pbc.pbc1(dy,box_len)
                dz = r[i][2] - r[j][2]
                dz = pbc.pbc1(dz,box_len)

                dist = (dx**2+dy**2+dz**2)**(1/2)    
            

                if dist < cutoff:
                
                    force[i][0] = force[i][0] + (48/(dist**14) - 24/(dist**6))*dx
                    force[i][1] = force[i][0] + (48/(dist**14) - 24/(dist**6))*dy
                    force[i][2] = force[i][0] + (48/(dist**14) - 24/(dist**6))*dz

                    force[j][0] = force[j][0] - (48/(dist**14) - 24/(dist**6))*dx
                    force[j][1] = force[j][1] - (48/(dist**14) - 24/(dist**6))*dy
                    force[j][2] = force[j][2] - (48/(dist**14) - 24/(dist**6))*dz

                    pot += 4*((1/dist)**12-(1/dist)**6) - 4*(1/cutoff**12 - 1/cutoff**6)
   
    return force, pot
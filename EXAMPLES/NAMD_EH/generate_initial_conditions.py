import numpy as np
import subprocess as sp

NAtoms = 6
NSamples = 500

Atom_labels  = []
GEOMS        = np.zeros(( NSamples, NAtoms, 3 ))
VELOCS       = np.zeros(( NSamples, NAtoms, 3 ))

INIT_FILE = open("initconds","r").readlines()
for count, line in enumerate(INIT_FILE):
    t = line.split()
    if ( len(t) == 2 and t[0] == "Index"):
        IND = int(t[1]) - 1
        offset = 2
        for at in range( NAtoms ):
            if (IND == 0):
                Atom_labels.append( INIT_FILE[count+offset+at].split()[0] )
            GEOMS[IND,at,:]     = np.array(INIT_FILE[count+offset+at].split()[2:5], dtype=float)
            VELOCS[IND,at,:]    = np.array(INIT_FILE[count+offset+at].split()[6:9], dtype=float)

# Make correct units G:[ANG] and V:[ANG / a.u.]
GEOMS  *= 0.529
VELOCS *= 1.000 # Not sure what these are yet

sp.call("rm -r TRAJ",shell=True)
sp.call("mkdir TRAJ",shell=True)

for traj in range( NSamples ):
    sp.call(f"mkdir TRAJ/traj-{traj}",shell=True)
    
    f1 = open(f"TRAJ/traj-{traj}/geometry_input.xyz","w")
    f2 = open(f"TRAJ/traj-{traj}/velocity_input.xyz","w")

    f1.write(f"{NAtoms}\nGeometry #{traj}\n")
    f2.write(f"{NAtoms}\nVelocity #{traj}\n")

    for at in range( NAtoms ):
        f1.write(  f"{Atom_labels[at]}  " + "  ".join(map("{:2.6f}".format, GEOMS[traj,at,:])) + "\n")
        f2.write(  f"{Atom_labels[at]}  " + "  ".join(map("{:2.6f}".format, VELOCS[traj,at,:])) + "\n")

    f1.close()
    f2.close()


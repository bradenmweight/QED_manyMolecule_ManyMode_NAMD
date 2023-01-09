import numpy as np
from random import gauss

import Eh

def get_random_force( DYN_PROPERTIES ):
            LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
            TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py
            NAtoms  = DYN_PROPERTIES["NAtoms"]
            masses = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ]) # TODO Change the shape in read_input.py
            dtI     = DYN_PROPERTIES["dtI"]
            #sigma = np.sqrt( 2 * LANGEVIN_LAMBDA * TEMP * masses ) # Gaussian width (NAtoms x 3)
            #xi    = np.array([gauss(0, 1) for j in range(3*NAtoms)]).reshape((NAtoms,3)) # (NAtoms x 3)
            #theta = np.array([gauss(0, 1) for j in range(3*NAtoms)]).reshape((NAtoms,3)) # (NAtoms x 3)
            #c     = 0.28867513459 # WHAT IS THIS ??? ~ BMW
            #xi    = np.array([gauss(0, sigma[j,d]/2) for j in range(NAtoms) for d in range(3)]).reshape((NAtoms,3)) # (NAtoms x 3)
            #theta = np.array([gauss(0, sigma[j,d]*c) for j in range(NAtoms) for d in range(3)]).reshape((NAtoms,3)) # (NAtoms x 3)

            #F_RAND = sigma * dtI**(3/2) * (xi/2 + c*theta)
            #F_RAND = sigma * dtI**(1/2) * (xi + theta) # dtI**2 outside function call
            

            
            return F_RAND

def get_damping_force( DYN_PROPERTIES ):
            LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
            masses = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ])
            F_DAMP          = -1.0 * LANGEVIN_LAMBDA * DYN_PROPERTIES["Atom_velocs_new"] * masses
            return F_DAMP

def Nuclear_X_Step(DYN_PROPERTIES):

    V       = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"] * 1.0

    # Propagate nuclear coordinates
    DYN_PROPERTIES["FORCE_NEW"] = Eh.get_Force(DYN_PROPERTIES) * 0.0
    a = DYN_PROPERTIES["FORCE_NEW"] / masses

    DYN_PROPERTIES["Atom_coords_new"] += V[:,:] * dtI + 0.5000000 * a[:,:] * dtI*dtI

    # FUNCTIONALITY FOR LANGEVIN DYNAMICS
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
            F_DAMP = get_damping_force( DYN_PROPERTIES )
            F_RAND = get_random_force( DYN_PROPERTIES )
            DYN_PROPERTIES["Atom_coords_new"] += (F_DAMP + F_RAND) / masses *  dtI*dtI

    return DYN_PROPERTIES

def Nuclear_V_Step(DYN_PROPERTIES):

    V       = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_velocs_old"] = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    DYN_PROPERTIES["FORCE_OLD"] = DYN_PROPERTIES["FORCE_NEW"] * 1.0 # Store old force
    DYN_PROPERTIES["Atom_velocs_old"] += DYN_PROPERTIES["Atom_velocs_new"]

    # Get new force
    DYN_PROPERTIES["FORCE_NEW"] = Eh.get_Force(DYN_PROPERTIES) * 0.0

    # Compute accelerations
    anew = DYN_PROPERTIES["FORCE_NEW"] / masses 
    aold = DYN_PROPERTIES["FORCE_OLD"] / masses

    DYN_PROPERTIES["Atom_velocs_new"] += 0.5000000 * (aold[:,:] + anew[:,:]) * dtI

    # FUNCTIONALITY FOR LANGEVIN DYNAMICS
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
            F_DAMP = get_damping_force( DYN_PROPERTIES )
            F_RAND = get_random_force( DYN_PROPERTIES )
            DYN_PROPERTIES["Atom_velocs_new"] += (F_DAMP + F_RAND) / masses *  dtI
            #B     = -dtI * LANGEVIN_LAMBDA * DYN_PROPERTIES["Atom_velocs_old"]
            #B    += sigma * np.sqrt(dtI) * xi - A * LANGEVIN_LAMBDA
            #DYN_PROPERTIES["Atom_velocs_new"] += B

    return DYN_PROPERTIES
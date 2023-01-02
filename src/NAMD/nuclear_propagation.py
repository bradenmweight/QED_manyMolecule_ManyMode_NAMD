import numpy as np

import Eh

def Nuclear_X_Step(DYN_PROPERTIES):

    X       = DYN_PROPERTIES["Atom_coords_new"] * 1.0
    V       = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"] * 1.0

    # Propagate nuclear coordinates
    DYN_PROPERTIES["FORCE_NEW"] = Eh.get_Force(DYN_PROPERTIES)
    a = DYN_PROPERTIES["FORCE_NEW"] / masses
    DYN_PROPERTIES["Atom_coords_new"] += V[:,:] * dtI + 0.5000000 * a[:,:] * dtI*dtI

    return DYN_PROPERTIES

def Nuclear_V_Step(DYN_PROPERTIES):

    V       = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_velocs_old"] = DYN_PROPERTIES["Atom_velocs_new"] * 1.0

    # Propagate nuclear coordinates
    DYN_PROPERTIES["FORCE_OLD"] = DYN_PROPERTIES["FORCE_NEW"] * 1.0 # Store old force
    DYN_PROPERTIES["FORCE_NEW"] = Eh.get_Force(DYN_PROPERTIES)

    anew = DYN_PROPERTIES["FORCE_NEW"] / masses 
    aold = DYN_PROPERTIES["FORCE_OLD"] / masses

    DYN_PROPERTIES["Atom_velocs_old"] += DYN_PROPERTIES["Atom_velocs_new"]
    DYN_PROPERTIES["Atom_velocs_new"] += 0.5000000 * (aold[:,:] + anew[:,:]) * dtI
    
    return DYN_PROPERTIES
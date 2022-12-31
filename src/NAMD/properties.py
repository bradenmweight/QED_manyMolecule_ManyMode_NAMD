import numpy as np

import Eh

def get_density_matrix( DYN_PROPERTIES ):
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.get_density_matrix(DYN_PROPERTIES)
    else:
        print("NAMD_METHOD not recognized. Quitting.")
        exit()


def compute_KE(DYN_PROPERTIES):
    KE = 0.0
    for at in range( DYN_PROPERTIES["NAtoms"] ):
        KE += 0.50000000 * DYN_PROPERTIES["MASSES"][at] * np.linalg.norm(DYN_PROPERTIES["Atom_velocs_new"][at,:])**2
    DYN_PROPERTIES["KE"] = KE
    return DYN_PROPERTIES

def compute_PE(DYN_PROPERTIES):
    PE = 0.0
    RHO = get_density_matrix(DYN_PROPERTIES)
    for state in range( DYN_PROPERTIES["NStates"] ):
        PE += RHO[state,state].real * DYN_PROPERTIES["DIAG_ENERGIES_NEW"][state]
    DYN_PROPERTIES["PE"] = PE
    return DYN_PROPERTIES


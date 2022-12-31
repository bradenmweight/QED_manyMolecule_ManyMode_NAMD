import numpy as np
import sys
import subprocess as sp
from time import time

import read_input
import nuclear_propagation
import output
import rotation

import Eh # Add other NAMD methods here

sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/")
sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/WFN_OVERLAP/PYTHON/")

import G16_TD


# This code with be the main control code for the NAMD
# We will make calls to electronic structure and
#   wavefunction overlaps here, which are taken care
#   of elsewhere
# Additionally, the mixed-quantum classical or semi-classical
#   dynamics will be handled elsewhere


def propagage_Mapping(DYN_PROPERTIES):
    """
    Wrapper for electronic propagation
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        Eh.propagage_Mapping(DYN_PROPERTIES)


def main( ):
    DYN_PROPERTIES = read_input.read()
    DYN_PROPERTIES = read_input.initialize_MD_variables(DYN_PROPERTIES)

    # Initialize electronic DOFs
    DYN_PROPERTIES = Eh.initialize_mapping(DYN_PROPERTIES)

    # Perform first electronic structure calculation
        # Get diagonal energies and gradients
    DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)
    output.save_data(DYN_PROPERTIES)

    # Start main MD loop
    for step in range( DYN_PROPERTIES["NSteps"] ):
        #print(f"Working on step {step} of { DYN_PROPERTIES['NSteps'] }")

        # Propagate nuclear coordinates
        DYN_PROPERTIES = nuclear_propagation.Nuclear_X_Step(DYN_PROPERTIES)

        # Perform jth electronic structure calculation
            # Get diagonal energies and grad
            # Get overlap and NACT
        DYN_PROPERTIES["MD_STEP"] += 1 # This needs to be exactly here for technical reasons.
        T0 = time()
        DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)
        print( "Total QM took %2.2f s." % (time() - T0) )

        if ( DYN_PROPERTIES["NStates"] >= 2 ):
            # Propagate electronic DOFs (Interpolated Hamiltonian)
            T0 = time()
            DYN_PROPERTIES = Eh.propagage_Mapping(DYN_PROPERTIES)
            print( "Electronic propagation took %2.2f s." % (time() - T0) )

            # Rotate electronic DOFs from t0 basis to t1 basis
            DYN_PROPERTIES = Eh.rotate_Mapping(DYN_PROPERTIES)

        # Propagate nuclear momenta
        DYN_PROPERTIES = nuclear_propagation.Nuclear_V_Step(DYN_PROPERTIES)
        
        """
        # NOT TESTED YET
        DYN_PROPERTIES = rotation.shift_COM(DYN_PROPERTIES)
        DYN_PROPERTIES = rotation.remove_rotations(DYN_PROPERTIES)
        """

        output.save_data(DYN_PROPERTIES)

if ( __name__ == "__main__" ):
    main()
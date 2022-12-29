import numpy as np
import sys
import subprocess as sp

import read_input
import Eh
import nuclear_propagation
import output

sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/")
sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/WFN_OVERLAP/PYTHON/")

import G16_TD


# This code with be the main control code for the NAMD
# We will make calls to electronic structure and
#   wavefunction overlaps here, which are taken care
#   of elsewhere
# Additionally, the mixed-quantum classical or semi-classical
#   dynamics will be handled elsewhere


def main( ):
    DYN_PROPERTIES = read_input.read()
    DYN_PROPERTIES = read_input.initialize_MD_variables(DYN_PROPERTIES)

    # Initialize electronic DOFs
    Eh.initialize_mapping(DYN_PROPERTIES)

    # Perform first electronic structure calculation
        # Get diagonal energies and gradients
    DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)

    print("R(Li-H) =", 0.529 * np.linalg.norm(DYN_PROPERTIES["Atom_coords_new"][0,:] - DYN_PROPERTIES["Atom_coords_new"][1,:]) )

    output.save_data(DYN_PROPERTIES)

    # Start main MD loop
    for step in range( DYN_PROPERTIES["NSteps"] ):
        #print(f"Working on step {step} of { DYN_PROPERTIES['NSteps'] }")

        # Propagate nuclear coordinates
        DYN_PROPERTIES = nuclear_propagation.Nuclear_X_Step(DYN_PROPERTIES)
        print("R(Li-H) =", 0.529 * np.linalg.norm(DYN_PROPERTIES["Atom_coords_new"][0,:] - DYN_PROPERTIES["Atom_coords_new"][1,:]) )
        #print("F(Li-H) =", DYN_PROPERTIES["FORCE_NEW"][1,:] - DYN_PROPERTIES["FORCE_NEW"][1,:] )


        # Perform jth electronic structure calculation
            # Get diagonal energies and grad
            # Get overlap and NACT
        DYN_PROPERTIES["MD_STEP"] +=1
        DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)

        # Propagate electronic DOFs (Interpolated Hamiltonian)

        # Transform electronic DOFs from t0 basis to t1 basis

        # Propagate nuclear momenta
        DYN_PROPERTIES = nuclear_propagation.Nuclear_V_Step(DYN_PROPERTIES)

        output.save_data(DYN_PROPERTIES)

if ( __name__ == "__main__" ):
    main()
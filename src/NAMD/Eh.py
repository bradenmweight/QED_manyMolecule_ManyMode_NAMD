import numpy as np

def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    z = np.zeros(( NStates ), dtype=complex)
    z[ISTATE] = 1.0 + 0.0j # Ehrenfest has no electronic sampling

    #DYN_PROPERTIES["ZPE"] = 0.0 # Ehrenfest has no ZPE 
    DYN_PROPERTIES["MAPPING_VARS"] = z

    return DYN_PROPERTIES

def get_Force(DYN_PROPERTIES):
    
    dEad    = DYN_PROPERTIES["DIAG_GRADIENTS"] # NStates x NAtoms x 3 (a.u.)
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    NStates = DYN_PROPERTIES["NStates"]

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        NACR    = DYN_PROPERTIES["NACR_APPROX_NEW"] # NStates x NStates x NAtoms x 3 (a.u.)
    else:
        NACR    = np.zeros(( NStates, NStates, NAtoms, 3 )) # NStates x NStates x NAtoms x 3 (a.u.)

    z = DYN_PROPERTIES["MAPPING_VARS"]

    # Do we need to find the state-independent force ?
    # It will already be included in the dEad terms at this point...
    #F0 = np.einsum("jad->ad", dEad[:,:,:]) / NStates

    F = np.zeros(( NAtoms, 3 ))
    rho = np.real( np.outer( np.conjugate(z), z ) )
    for j in range( NStates ):
        for k in range( j, NStates ):
            if ( j == k ):
                F[:,:] -= dEad[j,:,:] * rho[j,j]
                #print( dEad[j,-1,-1], rho[j,j] )
            else:
                F[:,:] -= -2 * NACR[j,k,:,:] * rho[j,k] # Double count upper triangle
    
    return F


if ( __name__ == "__main__" ):
    main()
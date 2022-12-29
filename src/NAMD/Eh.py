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

def propagage_Mapping(DYN_PROPERTIES):

    """
    Updates mapping variables with Verlet-like scheme
    TODO Add option for 4th-order Runge-Kutta
    """
    NStates = DYN_PROPERTIES["NStates"]
    z     = DYN_PROPERTIES["MAPPING_VARS"]

    Zreal = np.real(z)
    Zimag = np.imag(z)

    OVERLAP  = (DYN_PROPERTIES["OVERLAP_NEW"])[:-1,:-1] # Recall, we perform TD-DFT with one additional state

    print( "Wavefunction Overlap:" )
    print( np.round(OVERLAP,4) )

    Hamt0 = np.zeros(( NStates, NStates )) # t0 basis
    Hamt1 = np.zeros(( NStates, NStates )) # t1 basis

    #### t0 Ham ####
    # ADD NACT = dR/dt.NACR TO OFF-DIAGONALS
    if ( DYN_PROPERTIES["MD_STEP"] >= 2 ): 
        NACR_OLD = DYN_PROPERTIES["NACR_APPROX_OLD"]
        VELOC_OLD = DYN_PROPERTIES["Atom_velocs_old"]
    Ead_old   = DYN_PROPERTIES["DIAG_ENERGIES_OLD"]
    E_GS_t0   = Ead_old[0] * 1.0
    for j in range(NStates):
        for k in range(NStates):
            if ( j == k ):
                Hamt0[j,j] = Ead_old[j] - E_GS_t0
            else:
                if ( DYN_PROPERTIES["MD_STEP"] >= 2 ):
                    Ejk = Ead_old[j] - Ead_old[k]
                    Hamt0[j,k] = np.einsum("ad,ad->",NACR_OLD[j,k,:,:], VELOC_OLD[:,:] ) * Ejk

    #### t1 Ham in t0 basis ####
    NACR_NEW  = DYN_PROPERTIES["NACR_APPROX_NEW"]
    VELOC_NEW = DYN_PROPERTIES["Atom_velocs_new"]
    Ead_new   = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    for j in range(NStates):
        for k in range(NStates):
            if ( j == k ):
                Hamt1[j,j] = Ead_new[j] - E_GS_t0
            else:
                Ejk = Ead_new[j] - Ead_new[k]
                Hamt1[j,k] = np.einsum("ad,ad->",NACR_NEW[j,k,:,:], VELOC_NEW[:,:] ) * Ejk
        
    # TODO Check the direction of this rotation
    Hamt1 = OVERLAP.T @ Hamt1 @ OVERLAP # Rotate to t0 basis

    print(Hamt0)
    print(Hamt1)

    dtE    = DYN_PROPERTIES["dtE"]
    ESTEPS = DYN_PROPERTIES["ESTEPS"]

    for step in range( ESTEPS ):

        # Linear interpolation betwen t0 and t1
        H = Hamt0  + (step)/(ESTEPS) * ( Hamt1 - Hamt0 )

        # Propagate Imaginary first by dt/2
        Zimag -= 0.5000000 * H @ Zreal * dtE

        # Propagate Real by full dt
        Zreal += H @ Zimag * dtE
        
        # Propagate Imaginary final by dt/2
        Zimag -= 0.5000000 * H @ Zreal * dtE

    print( "Z (t0)", np.abs(DYN_PROPERTIES["MAPPING_VARS"]) )
    print( "Z (t1)", np.abs(Zreal + 1j*Zimag) )

    DYN_PROPERTIES["MAPPING_VARS"] = Zreal + 1j*Zimag

    return DYN_PROPERTIES

def rotate_Mapping(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    OVERLAP = (DYN_PROPERTIES["OVERLAP_NEW"])[:-1,:-1] # Recall, we perform TD-DFT with one additional state

    z = OVERLAP @ z # TODO Check the direction of this rotation
    #z = OVERLAP.T @ z # Maybe is this way

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    normalize_Mapping(DYN_PROPERTIES)

    return DYN_PROPERTIES


def normalize_Mapping(DYN_PROPERTIES): # Be careful with this. Ehrenfest it is okay.
    z = DYN_PROPERTIES["MAPPING_VARS"]
    norm = np.sum( np.abs(z) )
    print(f"Electronic Norm.: {np.round(norm,4)}")
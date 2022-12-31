import numpy as np
import random

import properties


def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]
    #DYN_PROPERTIES["ZPE"] = 0.0 # Ehrenfest has no ZPE 

    """
    ### MMST Style Initialization ###
    #z = np.zeros(( NStates ), dtype=complex)
    #z[ISTATE] = 1.0 + 0.0j # Ehrenfest has no electronic sampling
    """

    ### Spin-mapping Style Initialization ###
    Rw = 2*np.sqrt(NStates+1) # Radius of W Sphere

    # Initialize mapping radii
    r = np.zeros(( NStates )) # np.ones(( NStates )) * np.sqrt(gw)
    r[ISTATE] = np.sqrt( 2 )

    # Set real mapping variables
    z = np.zeros(( NStates ),dtype=complex)
    for i in range(NStates):
        phi = random.random() * 2 * np.pi # Azimuthal Angle -- Always Random
        z[i] = r[i] * ( np.cos( phi ) + 1j * np.sin( phi ) )
  
    # Check initial density matrix
    #rho = np.zeros((NStates,NStates),dtype=complex)
    #rho = 0.5 * np.outer( np.conjugate(z), z )
    #print("Initial Density Matrix:")
    #print( rho )

    DYN_PROPERTIES["MAPPING_VARS"] = z

    return DYN_PROPERTIES

def get_Force(DYN_PROPERTIES):
    
    dEad    = DYN_PROPERTIES["DIAG_GRADIENTS"] # NStates x NAtoms x 3 (a.u.)
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    NStates = DYN_PROPERTIES["NStates"]
    Ead     = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        NACR    = DYN_PROPERTIES["NACR_APPROX_NEW"] # NStates x NStates x NAtoms x 3 (a.u.)
    else:
        NACR    = np.zeros(( NStates, NStates, NAtoms, 3 )) # NStates x NStates x NAtoms x 3 (a.u.)

    z = DYN_PROPERTIES["MAPPING_VARS"]

    # Do we need to find the state-independent force ?
    # It will already be included in the dEad terms at this point...
    #F0 = np.einsum("jad->ad", dEad[:,:,:]) / NStates


    F = np.zeros(( NAtoms, 3 ))
    rho = np.real( properties.get_density_matrix(DYN_PROPERTIES) )
    for j in range( NStates ):
        for k in range( j, NStates ):
            if ( j == k ):
                F[:,:] -= 0.5 * dEad[j,:,:] * rho[j,j]
            else:
                Ejk = Ead[j] - Ead[k]
                F[:,:] -= 2 * 0.5 * rho[j,k] * NACR[j,k,:,:] * Ejk  # Double count upper triangle


    # FOR DEBUGGING CALCULATIONS, RHO = [1/2], equal super-position of ground and excited
    #F = np.zeros(( NAtoms, 3 ))
    #F[:,:] -= 0.5 * ( 0.5 * dEad[0,:,:] + 0.5 * dEad[1,:,:] ) + \
    #          0.5 * ( 0.5 * NACR[0,1,:,:] * (Ead[1] - Ead[0]) + 0.5 * NACR[1,0,:,:] * (Ead[0] - Ead[1]) )

    #print("FORCE:")
    #print(F[:,:])
    #print("NACR:")
    #print( NACR[0,1,:,:] )
    #print("dE")
    #print( (Ead[1] - Ead[0]) )
    #print("NACRx dE:")
    #print( NACR[0,1,:,:] * (Ead[0] - Ead[1]) )
    return F


def rotate_t0_to_t1(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap.
    if ( len(A.shape) == 1 ):
        #return S.T @ A
        return S @ A
    elif( len(A.shape) == 2 ):
        #return S.T @ A @ S
        return S @ A @ S.T
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numo array. Received {len(A.shape)}D array.")

def rotate_t1_to_t0(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap.
    if ( len(A.shape) == 1 ):
        #return S @ A
        return S.T @ A
    elif( len(A.shape) == 2 ):
        #return S @ A @ S.T
        return S.T @ A @ S
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numo array. Received {len(A.shape)}D array.")


def propagage_Mapping(DYN_PROPERTIES):
    NStates = DYN_PROPERTIES["NStates"]
    z       = DYN_PROPERTIES["MAPPING_VARS"]

    Zreal = np.real(z) * 1.0
    Zimag = np.imag(z) * 1.0

    OVERLAP  = (DYN_PROPERTIES["OVERLAP_NEW"]) # Recall, we perform TD-DFT with one additional state. Already removed.

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
            #else:
            #    if ( DYN_PROPERTIES["MD_STEP"] >= 2 ):
            #        Ejk = Ead_old[j] - Ead_old[k]
            #        Hamt0[j,k] = np.einsum("ad,ad->",NACR_OLD[j,k,:,:], VELOC_OLD[:,:] ) * Ejk

    #### t1 Ham in t0 basis ####
    NACR_NEW  = DYN_PROPERTIES["NACR_APPROX_NEW"]
    VELOC_NEW = DYN_PROPERTIES["Atom_velocs_new"]
    Ead_new   = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    for j in range(NStates):
        for k in range(NStates):
            if ( j == k ):
                Hamt1[j,j] = Ead_new[j] #- E_GS_t0 # Is it okay to subtract this here ? Identity will also rotate, no ? Maybe bad.
            #else:
            #    Ejk = Ead_new[j] - Ead_new[k]
            #    Hamt1[j,k] = np.einsum("ad,ad->",NACR_NEW[j,k,:,:], VELOC_NEW[:,:] ) * Ejk
        
    # TODO Check the direction of this rotation
    Hamt1 = rotate_t1_to_t0( DYN_PROPERTIES["OVERLAP_NEW"] , Hamt1 ) # Rotate to t0 basis

    Hamt1 -= np.identity(NStates) * E_GS_t0 # Subtract t0 reference 'after' rotation to t0 basis

    dtE    = DYN_PROPERTIES["dtE"]
    ESTEPS = DYN_PROPERTIES["ESTEPS"]

    if ( DYN_PROPERTIES["EL_PROP"] == "VV" ):
        """
        Propagate with second-order symplectic (Velocity-Verlet-like)
        """

        for step in range( ESTEPS ):

            # Linear interpolation betwen t0 and t1
            H = Hamt0 + (step)/(ESTEPS) * ( Hamt1 - Hamt0 )

            # Propagate Imaginary first by dt/2
            Zimag -= 0.5000000 * H @ Zreal * dtE

            # Propagate Real by full dt
            Zreal += H @ Zimag * dtE
            
            # Propagate Imaginary final by dt/2
            Zimag -= 0.5000000 * H @ Zreal * dtE

    elif ( DYN_PROPERTIES["EL_PROP"] == "RK" ):
        """
        Propagate with explicit 4th-order Runge-Kutta
        """

        def get_H( step, dt ):
            # Linear interpolation betwen t0 and t1
            H = Hamt0 + (step + dt/dtE)/(ESTEPS) * ( Hamt1 - Hamt0 )
            return H

        def f( y, H ):
            return -1j * H @ y

        yt = z.copy()

        for step in range( ESTEPS ):

            k1 = f(yt, get_H( step, 0 ))
            k2 = f(yt + k1*dtE/2, get_H( step, k1*dtE/2 ))
            k3 = f(yt + k2*dtE/2, get_H( step, k2*dtE/2 ))
            k4 = f(yt + k3*dtE, get_H( step, k3*dtE ))

            yt += 1/6 * ( k1 + 2*k2 + 2*k3 + k4 ) * dtE

        z = yt

    DYN_PROPERTIES["MAPPING_VARS"] = z

    return DYN_PROPERTIES

def rotate_Mapping(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    S = DYN_PROPERTIES["OVERLAP_NEW"]

    z = rotate_t0_to_t1( S, z )

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    DYN_PROPERTIES = normalize_Mapping(DYN_PROPERTIES)

    return DYN_PROPERTIES

def get_density_matrix( DYN_PROPERTIES ):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    return 0.500000 * np.outer( np.conjugate(z), z )


def normalize_Mapping(DYN_PROPERTIES): # Be careful with this. Ehrenfest it is okay I guess.
    z = DYN_PROPERTIES["MAPPING_VARS"]
    POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(DYN_PROPERTIES["NStates"])])
    norm = np.sum( POP )
    print(f"Electronic Norm.: {np.round(norm,4)} {np.round(norm,4)}")
    if ( abs(1.0 - norm) > 1e-5 and abs(1.0 - norm) < 1e-2 ):
        print(f"Electronic Norm.: {np.round(norm,4)}")
        print("Mapping norm is wrong.")
        #print("Renormalizing population.")
    elif( abs(1.0 - norm) > 1e-2 ):
        print(f"Electronic Norm.: {np.round(norm,4)}")
        print("Mapping norm is VERY wrong.")
        #print("Forcing renormalization... Should not need this...PLEASE FIX")
        #DYN_PROPERTIES["MAPPING_VARS"] /= norm1
        #exit()
    
    return DYN_PROPERTIES
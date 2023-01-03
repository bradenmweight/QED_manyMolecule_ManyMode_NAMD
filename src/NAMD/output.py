import numpy as np
import subprocess as sp

import properties

# MOVE THIS FUNCTION TO NEW FILE OUTPUT.py
def save_data(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]

    if ( DYN_PROPERTIES["MD_STEP"] == 0 ): 
        sp.call("rm -rf MD_OUTPUT ",shell=True)
        sp.call("mkdir MD_OUTPUT ",shell=True)

    with open("MD_OUTPUT/trajectory.xyz","a") as file01:
        file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
        file01.write(f"MD Step {DYN_PROPERTIES['MD_STEP']}\n")
        Atom_labels = DYN_PROPERTIES["Atom_labels"]
        Atom_coords = DYN_PROPERTIES["Atom_coords_new"] 
        for count, atom in enumerate( Atom_labels ):
            #file01.write(f"{atom}\t{Atom_coords[count,0]*0.529}\t{Atom_coords[count,1]*0.529}\t{Atom_coords[count,2]*0.529}\n")
            file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_coords[count,:]*0.529))  + "\n")

    with open("MD_OUTPUT/PES.dat","a") as file01:
        if ( NStates >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["DIAG_ENERGIES_NEW"]*27.2114 )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['DIAG_ENERGIES_NEW']*27.2114}\n" )

    with open("MD_OUTPUT/mapping_re.dat","a") as file01:
        if ( NStates >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].real )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['MAPPING_VARS'].real}\n" )

    with open("MD_OUTPUT/mapping_im.dat","a") as file01:
        if ( NStates >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].imag )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['MAPPING_VARS'].imag}\n" )

    with open("MD_OUTPUT/Population.dat","a") as file01:
        if ( NStates >= 2 ):
            POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(NStates)])
            PSUM = np.sum(POP)
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,POP )) + "  %2.8f" % (PSUM) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  1.0 1.0\n" )

    if ( NStates >= 2 ):
        with open("MD_OUTPUT/Coherence_re.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
            RHO = np.real(properties.get_density_matrix( DYN_PROPERTIES ))
            COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)]).real
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

        with open("MD_OUTPUT/Coherence_im.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
            RHO = np.real(properties.get_density_matrix( DYN_PROPERTIES ))
            COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)]).imag
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

    with open("MD_OUTPUT/Energy.dat","a") as file01:
        
        DYN_PROPERTIES = properties.compute_KE(DYN_PROPERTIES)
        DYN_PROPERTIES = properties.compute_PE(DYN_PROPERTIES)

        KE = DYN_PROPERTIES["KE"] * 27.2114
        PE = DYN_PROPERTIES["PE"] * 27.2114
        TE = KE + PE

        file01.write(f"{DYN_PROPERTIES['MD_STEP']}  " + "%2.6f  %2.6f  %2.6f\n" % (KE,PE,TE))


        with open("MD_OUTPUT/Overlap.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                OVERLAP = DYN_PROPERTIES['OVERLAP_NEW'] * 1.0
                OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

        with open("MD_OUTPUT/NACT.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                NACT = DYN_PROPERTIES['NACT_NEW'] * 1.0
                NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )

        with open("MD_OUTPUT/Overlap_uncorrected.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                OVERLAP = DYN_PROPERTIES['OVERLAP_NEW_uncorrected'] * 1.0
                OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

        with open("MD_OUTPUT/NACT_uncorrected.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                NACT = DYN_PROPERTIES['NACT_NEW_uncorrected'] * 1.0
                NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )
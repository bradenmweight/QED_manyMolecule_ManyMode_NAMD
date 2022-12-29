import numpy as np
import subprocess as sp

import properties

# MOVE THIS FUNCTION TO NEW FILE OUTPUT.py
def save_data(DYN_PROPERTIES):
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
        if ( DYN_PROPERTIES["NStates"] >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["DIAG_ENERGIES_NEW"]*27.2114 )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['DIAG_ENERGIES_NEW']*27.2114}\n" )

    with open("MD_OUTPUT/mapping_re.dat","a") as file01:
        if ( DYN_PROPERTIES["NStates"] >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].real )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['MAPPING_VARS'].real}\n" )

    with open("MD_OUTPUT/mapping_im.dat","a") as file01:
        if ( DYN_PROPERTIES["NStates"] >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].imag )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['MAPPING_VARS'].imag}\n" )

    with open("MD_OUTPUT/Energy.dat","a") as file01:
        
        DYN_PROPERTIES = properties.compute_KE(DYN_PROPERTIES)
        DYN_PROPERTIES = properties.compute_PE(DYN_PROPERTIES)

        KE = DYN_PROPERTIES["KE"] * 27.2114
        PE = DYN_PROPERTIES["PE"] * 27.2114
        TE = KE + PE

        file01.write(f"{DYN_PROPERTIES['MD_STEP']}  " + "%2.6f  %2.6f  %2.6f\n" % (KE,PE,TE))







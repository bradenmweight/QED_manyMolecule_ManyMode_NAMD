import numpy as np
import subprocess as sp


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
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["DIAG_ENERGIES"]*27.2114 )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {DYN_PROPERTIES['DIAG_ENERGIES']*27.2114}\n" )


    #with open("KE.dat","a") as file01:
    #    file01.write(f"{compute_KE(DYN_PROPERTIES)}\n")
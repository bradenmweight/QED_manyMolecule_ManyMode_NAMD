import numpy as np


def read():

    DYN_PROPERTIES = {} # Everything will come from here

    # Read input file
    input_lines = open('NAMD.in','r').readlines()
    for count, line in enumerate(input_lines):
        ### Clean line and Check for comments ###
        t = line.split()
        if ( len(t) == 0 or line.split()[0] in ["#","!"] ): continue # Check for comment line
        t = [ j.strip() for j in line.split("=") ]
        tnew = []
        for tj in t:
            if ( "#" in tj.split() ):
                tnew.append( tj.split()[:tj.split().index("#")][0] )
                break
            if ( "!" in tj.split() ):
                tnew.append( tj.split()[:tj.split().index("!")][0] )
                break
            else:
                tnew.append(tj)
        t = tnew

        if ( len(t) == 2 ):
            
            # Look for NStates
            if ( t[0].lower() == "nstates" ):
                try:
                    DYN_PROPERTIES["NStates"] = int( t[1] )
                except ValueError:
                    print(f"\t'NSteps' must be an integer: '{t[1]}'")
                    exit()

            # Look for NSteps
            if ( t[0].lower() == "nsteps" ):
                try:
                    DYN_PROPERTIES["NSteps"] = int( t[1] )
                except ValueError:
                    print(f"\t'NSteps' must be an integer: '{t[1]}'")
                    exit()

            # Look for ISTATE
            if ( t[0].lower() == "istate" ):
                try:
                    DYN_PROPERTIES["ISTATE"] = int( t[1] )
                except ValueError:
                    print(f"\t'ISTATE' must be an integer: '{t[1]}'")
                    exit()

            # Look for FUNCTIONAL
            if ( t[0].lower() == "FUNCTIONAL".lower() ):
                DYN_PROPERTIES["FUNCTIONAL"] = t[1]
                # The accuracy of input is up to the user...scary!

            # Look for BASIS
            if ( t[0].lower() == "BASIS".lower() ):
                DYN_PROPERTIES["BASIS_SET"] = t[1]
                # The accuracy of input is up to the user...scary!

            # Look for ISTATE
            if ( t[0].lower() == "dtI".lower() ):
                try:
                    DYN_PROPERTIES["dtI"] = float( t[1] ) * 41.341 # 41.341 a.u. / fs
                except ValueError:
                    print(f"\t'dtI' must be a float: '{t[1]}'")
                    exit()
            
            # Look for ESTEPS
            if ( t[0].lower() == "ESTEPS".lower() ):
                try:
                    DYN_PROPERTIES["ESTEPS"] = int( t[1] )
                except ValueError:
                    print(f"\t'ESTEPS' must be an integer: '{t[1]}'")
                    exit()

            # Look for NCPUS
            if ( t[0].lower() == "NCPUS".lower() ):
                try:
                    DYN_PROPERTIES["NCPUS"] = int( t[1] )
                except ValueError:
                    print(f"\t'NCPUS' must be an integer: '{t[1]}'")
                    exit()

            # Look for MEMORY
            if ( t[0].lower() == "MEMORY".lower() ):
                try:
                    DYN_PROPERTIES["MEMORY"] = int( t[1] )
                except ValueError:
                    print(f"\t'MEMORY' must be an integer: '{t[1]}'")
                    exit()

            # Look for CHARGE
            if ( t[0].lower() == "CHARGE".lower() ):
                try:
                    DYN_PROPERTIES["CHARGE"] = int( t[1] )
                except ValueError:
                    print(f"\t'CHARGE' must be an integer: '{t[1]}'")
                    exit()

            # Look for CHARGE
            if ( t[0].lower() == "MULTIPLICITY".lower() ):
                try:
                    DYN_PROPERTIES["MULTIPLICITY"] = int( t[1] )
                except ValueError:
                    print(f"\t'MULTIPLICITY' must be an integer: '{t[1]}'")
                    exit()

            # Look for VELOC
            if ( t[0].lower() == "VELOC".lower() ):
                DYN_PROPERTIES["VELOC"] = t[1]
                # Later, we will check this. We will automatically generate Wigner is VELOC != read

            # Look for RUN_ELEC_STRUC
            if ( t[0].lower() == "RUN_ELEC_STRUC".lower() ):
                DYN_PROPERTIES["RUN_ELEC_STRUC"] = t[1]
                # The accuracy of input is up to the user...scary!

            # Look for SBATCH_G16
            if ( t[0].lower() == "SBATCH_G16".lower() ):
                DYN_PROPERTIES["SBATCH_G16"] = t[1]
                # The accuracy of input is up to the user...scary!

        else:
            print( f"Error: Input is wrong at line {count+1}: {line}" )
            print( f"\tToo many '='" )
            exit()

    try:
        print("MANDATORY INPUT VARIABLES:")
        print( "\tNStates =", DYN_PROPERTIES["NStates"] )
        print( "\tNSteps =", DYN_PROPERTIES["NSteps"] )
        print( "\tdtI =", DYN_PROPERTIES["dtI"]/41.341, "(fs)" )
        print( "\tESTEPS =", DYN_PROPERTIES["ESTEPS"] )
        print( "\tISTATE =", DYN_PROPERTIES["ISTATE"] )
        print( "\tFUNCTIONAL =", DYN_PROPERTIES["FUNCTIONAL"] )
        print( "\tBASIS_SET =", DYN_PROPERTIES["BASIS_SET"] )
        print( "\tCHARGE =", DYN_PROPERTIES["CHARGE"] )
        print( "\tMULTIPLICITY =", DYN_PROPERTIES["MULTIPLICITY"] )
        print( "\tMEMORY =", DYN_PROPERTIES["MEMORY"], "(GB)" )
    except KeyError:
        print("Input file is missing mandatory entries. Check it.")

    assert( DYN_PROPERTIES["ISTATE"] <= DYN_PROPERTIES["NStates"]-1 ), "ISTATE must be less than the total number of states."


    return DYN_PROPERTIES



def read_geom():
    """
    TODO Add checks for XYZ user input
    """
    XYZ_File = open("geometry_input.xyz","r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_labels = []
    Atom_coords_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_labels.append( t[0] )
        Atom_coords_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ]) / 0.529 # Ang -> a.u.

    return Atom_labels, Atom_coords_new

def set_masses(Atom_labels):
    mass_amu_to_au = 1837/1.007 # au / amu
    masses_amu = { "H":1.007, 
                   "He":4.00260,
                   "Li":6.941,
                   "Be":9.01218,
                   "B":10.81,
                   "C":12.011,
                   "N":14.0067,
                   "O":15.9994,
                   "F":18.998403,
                   "Ne":20.179}

    masses = []
    for at in Atom_labels:
        masses.append( masses_amu[at] )
    return np.array(masses) * mass_amu_to_au

def get_initial_velocs(DYN_PROPERTIES):
    
    Atom_labels = DYN_PROPERTIES["Atom_labels"]
    masses      = DYN_PROPERTIES["MASSES"]

    # TODO Get Wigner distribution for initial velocities

    # For now, choose Maxwell-Boltzmann Distribution
    if (True):
        import random
        velocs = np.zeros(( len(Atom_labels), 3 ))
        T = 10 # K
        kT  = T * (0.025/300) / 27.2114 # K -> eV -> au
        V0  = np.sqrt( 2 * kT / masses ) # Average zero velocity
        SIG = kT / masses
        for at,atom in enumerate(Atom_labels):
            for d in range(3):
                velocs[at,d] = random.gauss( V0[at], SIG[at] )
    else:
        velocs = np.zeros(( len(Atom_labels), 3 ))

    return velocs

def initialize_MD_variables(DYN_PROPERTIES):
    
    DYN_PROPERTIES["Atom_labels"], DYN_PROPERTIES["Atom_coords_new"] = read_geom()
    DYN_PROPERTIES["MD_STEP"] = 0
    DYN_PROPERTIES["NAtoms"] = len( DYN_PROPERTIES["Atom_labels"] )

    DYN_PROPERTIES["MASSES"] = set_masses(DYN_PROPERTIES["Atom_labels"])

    DYN_PROPERTIES["Atom_velocs_new"] = get_initial_velocs(DYN_PROPERTIES)

    return DYN_PROPERTIES
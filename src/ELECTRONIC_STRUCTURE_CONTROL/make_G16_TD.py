import numpy as np
import subprocess as sp
import os, sys
import time

import get_cartesian_gradients
import get_diagonal_electronic_energies

sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/WFN_OVERLAP/PYTHON/")
import G16_NAC

# Obtain a geometry

# Move older directories: GS_NEW -> GS_OLD, TD_NEW -> TD_OLD

# Remove OVERLAP directory

# Make directories for GS_NEW (with force), TD_NEW (S1 force), OVERLAP

# Submit all three jobs

# Once TD is done, submit TD for all forces (Sj >= S2) simultaneously (OpenMP used here ?)

# Run data-collecting scripts in top directory (Need to add GS directories to this procedure yet.) # TODO

def check_geometry(Atom_labels,Atom_coords_new):

    assert ( isinstance(Atom_labels, list) ), "Atoms labels needs to be a list" 
    assert ( isinstance(Atom_labels[0], str) ), "Atoms labels needs to be a list of strings"
    assert ( isinstance(Atom_coords_new, type(np.array([])) ) ) , "Atoms coordinates need to be numpy array"

def clean_directory(NStates, MD_STEP):

    sp.call( "rm -rf GS_OLD" ,shell=True)
    sp.call( "rm -rf TD_OLD_S1" ,shell=True)
    if ( MD_STEP >= 2 ): sp.call( "rm -rf DIMER" ,shell=True)
    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            sp.call( f"rm -rf TD_OLD_S{state}" ,shell=True)

    if ( os.path.exists("TD_NEW_S1") ):
        sp.call( "mv TD_NEW_S1 TD_OLD_S1" ,shell=True)
        sp.call( "mv GS_NEW GS_OLD" ,shell=True)
        if ( NStates >= 3 ):
            for state in range( 2, NStates ):
                sp.call( f"mv TD_NEW_S{state} TD_OLD_S{state}" ,shell=True)

    sp.call( "mkdir -p GS_NEW" ,shell=True)
    sp.call( "mkdir -p TD_NEW_S1" ,shell=True)
    if ( MD_STEP >= 1 ): sp.call( "mkdir -p DIMER" ,shell=True)

    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            sp.call( f"mkdir -p TD_NEW_S{state}" ,shell=True)

    # ADD DIRECTORY TO KEEP TRACK OF PREVIOUS PHASE INFORMATION. ADD LATER. # TODO


def generate_inputs(DYN_PROPERTIES):

    NStates         = DYN_PROPERTIES["NStates"]
    Atom_labels     = DYN_PROPERTIES["Atom_labels"]
    Atom_coords_new = DYN_PROPERTIES["Atom_coords_new"]
    FUNCTIONAL      = DYN_PROPERTIES["FUNCTIONAL"]
    BASIS_SET       = DYN_PROPERTIES["BASIS_SET"]
    MEM             = DYN_PROPERTIES["MEMORY"]
    NCPUS           = DYN_PROPERTIES["NCPUS"]
    CHARGE          = DYN_PROPERTIES["CHARGE"]
    MULTIPLICITY    = DYN_PROPERTIES["MULTIPLICITY"]
    MD_STEP         = DYN_PROPERTIES["MD_STEP"]
    RUN_ELEC_STRUC  = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    SBATCH_G16      = DYN_PROPERTIES["SBATCH_G16"]

    #assert(FUNCTIONAL.upper() not in ['AM1', 'PM3', 'PM6', 'PM7'] ), "CI Overlap code does not work with semi-empirical Hamiltonians."

    def write_header(file01,MEM,NCPUS):
        file01.write(f"%chk=geometry.chk\n")
        file01.write(f"%nprocshared={NCPUS}\n")
        file01.write(f"%mem={MEM}GB\n\n")

    def write_geom(file01, Atom_labels, Atom_coords_new, MD_STEP, CHARGE, MULTIPLICITY, method=[None,None]):
        file01.write(f"MD Step {MD_STEP}\n\n")
        file01.write(f"{CHARGE} {MULTIPLICITY}\n")
        for at in range( len(Atom_labels) ):
            file01.write( "%s  %2.8f  %2.8f  %2.8f \n" % (Atom_labels[at],Atom_coords_new[at,0],Atom_coords_new[at,1],Atom_coords_new[at,2]) )
        if ( method[0] == "DIMER" ):
            coords_dimer = method[1]
            for at in range( len(Atom_labels) ):
                file01.write( "%s  %2.8f  %2.8f  %2.8f \n" % (Atom_labels[at],coords_dimer[at,0],coords_dimer[at,1],coords_dimer[at,2]) )
        file01.write("\n\n\n\n\n\n\n\n")

    # Ground State for New Geometry
    os.chdir("GS_NEW/")
    file01 = open("geometry.com","w")
    if ( MD_STEP > 1 ): 
        file01.write(f"%oldchk=../GS_OLD/geometry.chk\n")
    write_header(file01,MEM,NCPUS)
    file01.write(f"# {FUNCTIONAL}/{BASIS_SET} FORCE nosym pop=full\n\n") ### MAIN LINE ###
    write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
    os.chdir("../")

    # Excited State for New Geometry (Use converged wavefunctions from GS at same geometry)
    os.chdir("TD_NEW_S1/")
    file01 = open("geometry.com","w")
    file01.write(f"%oldchk=../GS_NEW/geometry.chk\n")
    write_header(file01,MEM,NCPUS)
    # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
    file01.write(f"# {FUNCTIONAL}/{BASIS_SET} TD=(singlets,nstates={NStates},root=1) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
    write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
    os.chdir("../")

    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            # Excited State for New Geometry (Use converged wavefunctions from TD root=1)
            os.chdir(f"TD_NEW_S{state}/")
            file01 = open("geometry.com","w")
            file01.write(f"%oldchk=../TD_NEW_S1/geometry.chk\n")
            write_header(file01,MEM,NCPUS)
            # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} TD=(read,singlets,nstates={NStates},root={state}) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
            write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
            os.chdir("../")
    
    # Dimer method for atomic orbital overlaps
    if ( MD_STEP >= 1 ):
        Atom_coords_old = DYN_PROPERTIES["Atom_coords_old"]
        os.chdir("DIMER/")
        file01 = open("geometry.com","w")
        write_header(file01,MEM,NCPUS)
        if ( FUNCTIONAL.upper() in ['AM1', 'PM3', 'PM6', 'PM7'] ): # By default, gaussian does not compute AO overlap for semi-empirical Hamiltonians
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} IOp(3/41=2000) nosymm iop(2/12=3,3/33=1) guess=only pop=full\n\n") ### MAIN LINE ###
        else:
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} nosymm iop(2/12=3,3/33=1) guess=only pop=full\n\n") ### MAIN LINE ###
        write_geom(file01,Atom_labels,Atom_coords_old,MD_STEP,CHARGE,MULTIPLICITY,method=["DIMER",Atom_coords_new])
        os.chdir("../")



def submit_jobs(DYN_PROPERTIES):
    NStates         = DYN_PROPERTIES["NStates"]
    RUN_ELEC_STRUC  = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    SBATCH_G16      = DYN_PROPERTIES["SBATCH_G16"]
    MD_STEP      = DYN_PROPERTIES["MD_STEP"]

    def submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP):
        if ( RUN_ELEC_STRUC == "SUBMIT_SBATCH" ):
            sp.call(f"cp {SBATCH_G16} .", shell=True)
            sp.call(f"sbatch {SBATCH_G16.split('/')[-1]}", shell=True)
            t0 = time.time()
            sleep_time = 0
            sleep_limit = 60 * 20 # 20 minutes, if electronic structure takes longer, we should not do NAMD
            sleep_check = 1 # Check gaussian output every {} seconds
            while ( True ):
                time.sleep(sleep_check)
                sleep_time += sleep_check # Add 5 seconds to sleep timer
                if ( sleep_limit > sleep_limit ):
                    print(f"\tWARNING! Gaussian did not finish normally in the following directory:\n{os.getcwd()}", )
                    exit()
                try:
                    check = open("geometry.out","r").readlines()[-1]
                except FileNotFoundError:
                    continue
                if ( check.split()[:4] == "Normal termination of Gaussian".split() ):
                    print("\tGaussian terminated normally in %2.2f s. (%s)" % (time.time() - t0, os.getcwd().split("/")[-1]) )
                    sp.call(f"formchk geometry.chk > /dev/null 2>&1", shell=True) # For dipole matrix calculations
                    break

        elif( RUN_ELEC_STRUC == "USE_CURRENT_NODE" ):
            t0 = time.time()
            sp.call(f"g16 < geometry.com > geometry.out", shell=True, stdout=True)
            check = open("geometry.out","r").readlines()[-1]
            if ( check.split()[:4] == "Normal termination of Gaussian".split() ):
                print("\tGaussian terminated normally in %2.2f s. (%s)" % (time.time() - t0, os.getcwd().split("/")[-1]) )
                sp.call(f"formchk geometry.chk > /dev/null 2>&1", shell=True) # For dipole matrix calculations
            else:
                print(f"\tWARNING! Gaussian did not finish normally in the following directory:\n{os.getcwd()}", )
        else:
            # This is a test. Do nothing 
            print(f"Testing. I will not submit/run electronic structure calculations for step {MD_STEP}.")


    print(f"Submitting electronic structure for step {MD_STEP}.")

    # GS and first ES must to be serial jobs
    os.chdir("GS_NEW/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
    os.chdir("../")

    os.chdir("TD_NEW_S1/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
    os.chdir("../")

    ### ADD PARALLELIZATION HERE FOR COMPUTING ALL THE EXCITED STATE FORCES IF NSTATES >= 3 ###

    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            os.chdir(f"TD_NEW_S{state}/")
            submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
            os.chdir("../")

    os.chdir("DIMER/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
    os.chdir("../")

def main(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"] # Total number of electronic states
    Atom_labels = DYN_PROPERTIES["Atom_labels"]
    Atom_coords_new = DYN_PROPERTIES["Atom_coords_new"]
    MD_STEP = DYN_PROPERTIES["MD_STEP"]
    
    if ( not os.path.exists("G16") ):
        sp.call("mkdir G16", shell=True)
    os.chdir("G16")
    
    check_geometry(Atom_labels,Atom_coords_new)
    clean_directory(NStates,MD_STEP)
    generate_inputs(DYN_PROPERTIES)
    submit_jobs(DYN_PROPERTIES)

    os.chdir("../")

    get_cartesian_gradients.main(DYN_PROPERTIES)
    get_diagonal_electronic_energies.main(DYN_PROPERTIES)
    if ( MD_STEP >= 1 ):
        G16_NAC.main(DYN_PROPERTIES)


def read_XYZ():
    XYZ_File = open("geometry_new.xyz","r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_labels = []
    Atom_coords_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_labels.append( t[0] )
        Atom_coords_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ])

    return Atom_labels, Atom_coords_new

if ( __name__ == "__main__" ):

    # THIS CODE NEEDS TO BE RUN TWICE WHEN USING THE EXAMPLE

    sp.call("module load gaussian", shell=True)
    sp.call("module load intel", shell=True)

    Atom_labels, Atom_coords_new = read_XYZ()

    DYN_PROPERTIES = {"Atom_labels":Atom_labels, "Atom_coords_new":Atom_coords_new }
    Atom_coords_old = Atom_coords_new * 1.0
    Atom_coords_old[0,0] += 0.005
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"]

    DYN_PROPERTIES["NStates"] = 4
    DYN_PROPERTIES["NAtoms"] = len(DYN_PROPERTIES["Atom_labels"])
    DYN_PROPERTIES["FUNCTIONAL"] = "AM1"
    DYN_PROPERTIES["CHARGE"] = 0
    DYN_PROPERTIES["MULTIPLICITY"] = 1
    DYN_PROPERTIES["BASIS_SET"] = "sto-3g"
    DYN_PROPERTIES["MEMORY"] = 5
    DYN_PROPERTIES["NCPUS"] = 1
    DYN_PROPERTIES["MD_STEP"] = 1
    DYN_PROPERTIES["RUN_ELEC_STRUC"] = "SUBMIT_SBATCH" # "SUBMIT_SBATCH", "USE_CURRENT_NODE", "TEST"
    DYN_PROPERTIES["SBATCH_G16"] = "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_INPUT/EXAMPLE/submit.gaussian" # For "SUBMIT_SBATCH" in previous only
    main(DYN_PROPERTIES)

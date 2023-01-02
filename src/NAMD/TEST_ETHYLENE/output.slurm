Traceback (most recent call last):
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src//NAMD/NAMD_main.py", line 86, in <module>
    main()
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src//NAMD/NAMD_main.py", line 59, in main
    DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/G16_TD.py", line 277, in main
    DYN_PROPERTIES = get_cartesian_gradients.main(DYN_PROPERTIES)
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/get_cartesian_gradients.py", line 52, in main
    read_Gradients(DIAG_GRADIENTS,DYN_PROPERTIES)
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/get_cartesian_gradients.py", line 35, in read_Gradients
    DIAG_GRADIENTS[1,:,:] = read_FCHK_GRAD( DYN_PROPERTIES["NAtoms"] )
  File "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/get_cartesian_gradients.py", line 12, in read_FCHK_GRAD
    lines = open("geometry.fchk","r").readlines()
FileNotFoundError: [Errno 2] No such file or directory: 'geometry.fchk'

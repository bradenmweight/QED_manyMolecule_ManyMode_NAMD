%oldchk=../TD_NEW_S1/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/sto-3g TD=(read,singlets,nstates=4,root=2) FORCE nosym pop=full guess=read

MD Step 1

0 1
Li  0.00000000  0.00000000  0.00000000 
H  0.00000000  0.00000000  1.00000000 









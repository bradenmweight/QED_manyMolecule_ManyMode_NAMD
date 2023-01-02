%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/sto-3g TD=(singlets,nstates=4,root=1) FORCE nosym pop=full guess=read

MD Step 1

0 1
Li  0.00000000  0.00000000  0.00000000 
H  0.00000000  0.00000000  1.00000000 









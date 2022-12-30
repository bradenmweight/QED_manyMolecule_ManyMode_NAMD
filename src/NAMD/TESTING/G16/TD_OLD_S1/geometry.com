%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/3-21G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 0

0 1
Li  0.00000000  0.00000000  0.00000000 
F  0.00000000  0.00000000  1.00000000 









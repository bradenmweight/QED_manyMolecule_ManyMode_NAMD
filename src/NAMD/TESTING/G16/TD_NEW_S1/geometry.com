%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/3-21G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1

0 1
Li  0.00833099  0.00833316  -0.02813557 
F  0.00503465  0.00503488  1.01835926 









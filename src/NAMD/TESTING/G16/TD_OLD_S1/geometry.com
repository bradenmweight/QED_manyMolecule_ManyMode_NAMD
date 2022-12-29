%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 2

0 1
H  0.02278707  0.68917313  -0.65651738 
H  0.02286242  0.69488250  0.70774618 
O  0.00536492  -0.01593397  0.00501499 









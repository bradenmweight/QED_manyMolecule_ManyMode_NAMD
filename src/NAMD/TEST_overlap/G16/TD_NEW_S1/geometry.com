%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 0

0 1
C  -0.38329520  0.73226544  0.00000000 
H  0.14986855  -0.19543948  0.00000000 
H  -1.45329520  0.73226544  0.00000000 
C  0.29197911  1.90724273  0.00000000 
H  -0.24118464  2.83494765  0.00000000 
H  1.36197911  1.90724273  0.00000000 









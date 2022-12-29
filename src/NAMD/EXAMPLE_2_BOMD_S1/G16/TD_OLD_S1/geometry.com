%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 4

0 1
Li  0.00610353  0.00610291  0.02436219 
H  0.01584225  0.01584182  2.88998647 









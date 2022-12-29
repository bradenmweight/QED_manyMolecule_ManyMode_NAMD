%oldchk=../TD_NEW_S1/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/STO-3G SCF=XQC TD=(read,singlets,nstates=4,root=2) FORCE nosym pop=full guess=read

MD Step 1000

0 1
Li  1.76309086  1.70981190  2.34273130 
H  2.32428623  2.69174964  1.32981945 









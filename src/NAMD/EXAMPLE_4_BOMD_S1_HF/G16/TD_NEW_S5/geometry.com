%oldchk=../TD_NEW_S1/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/3-21G* SCF=XQC TD=(read,singlets,nstates=8,root=5) FORCE nosym pop=full guess=read

MD Step 1000

0 1
Li  2.21054800  2.20729076  1.69947539 
F  0.66750128  0.66861839  3.85415908 









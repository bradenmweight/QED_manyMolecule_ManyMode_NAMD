%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/3-21G* SCF=XQC FORCE nosym pop=full

MD Step 1000

0 1
Li  2.21054800  2.20729076  1.69947539 
F  0.66750128  0.66861839  3.85415908 









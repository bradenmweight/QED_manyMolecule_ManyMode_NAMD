%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 1000

0 1
Li  1.82924668  1.82934944  2.56831455 
H  1.86855023  1.86876606  -0.22475364 









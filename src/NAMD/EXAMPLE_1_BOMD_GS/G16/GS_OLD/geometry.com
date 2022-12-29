%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 999

0 1
Li  1.82692732  1.82702943  2.57701163 
H  1.87005992  1.87027934  -0.29917846 









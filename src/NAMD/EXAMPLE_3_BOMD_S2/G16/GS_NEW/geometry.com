%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 1000

0 1
Li  1.76309086  1.70981190  2.34273130 
H  2.32428623  2.69174964  1.32981945 









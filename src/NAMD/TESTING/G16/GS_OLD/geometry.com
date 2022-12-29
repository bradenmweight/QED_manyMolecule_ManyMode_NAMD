%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 10

0 1
Li  0.00919629  0.00919615  0.05058556 
H  0.02347538  0.02347770  2.73818950 









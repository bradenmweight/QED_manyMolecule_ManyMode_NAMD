%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 999

0 1
Li  2.05784675  1.82840903  2.28487710 
H  0.27938737  1.85966853  1.71367059 








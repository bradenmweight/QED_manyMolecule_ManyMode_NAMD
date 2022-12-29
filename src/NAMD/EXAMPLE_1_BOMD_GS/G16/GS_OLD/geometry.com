%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 32

0 1
Li  0.03384584  0.03384126  0.06913311 
H  0.22997601  0.23003715  2.98677867 









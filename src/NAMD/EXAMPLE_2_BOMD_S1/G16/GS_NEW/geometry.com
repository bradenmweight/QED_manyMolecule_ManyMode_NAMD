%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 5

0 1
Li  0.00764420  0.00764342  0.03627239 
H  0.01970090  0.01970037  2.82236966 









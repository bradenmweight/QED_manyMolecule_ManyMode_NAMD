%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC FORCE nosym pop=full

MD Step 3

0 1
H  0.03507766  0.83170462  -0.78873981 
H  0.03563322  0.85951591  0.88639707 
O  0.00791387  -0.03108617  0.00622062 









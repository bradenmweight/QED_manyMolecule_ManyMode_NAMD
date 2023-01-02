%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/sto-3g nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 1

0 1
Li  0.00000000  0.00000000  0.00000000 
H  0.00000000  0.00000000  1.00000000 
Li  0.00000000  0.00000000  0.00000000 
H  0.00000000  0.00000000  1.00000000 









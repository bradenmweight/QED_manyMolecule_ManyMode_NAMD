%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/sto-3g FORCE nosym pop=full

MD Step 1

0 1
Li  0.00000000  0.00000000  0.00000000 
H  0.00000000  0.00000000  1.00000000 









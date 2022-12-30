%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/3-21G nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 1

0 1
Li  0.00000000  0.00000000  0.00000000 
F  0.00000000  0.00000000  1.00000000 
Li  0.00833099  0.00833316  -0.02813557 
F  0.00503465  0.00503488  1.01835926 









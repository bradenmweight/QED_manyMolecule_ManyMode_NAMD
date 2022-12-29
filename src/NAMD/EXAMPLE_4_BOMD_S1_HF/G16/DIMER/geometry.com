%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/3-21G* nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 1000

0 1
Li  2.21337501  2.21008278  1.69066196 
F  0.66499332  0.66612330  3.85590398 
Li  2.21054800  2.20729076  1.69947539 
F  0.66750128  0.66861839  3.85415908 









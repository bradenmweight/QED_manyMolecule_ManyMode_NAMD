%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 5

0 1
Li  0.00610353  0.00610291  0.02436219 
H  0.01584225  0.01584182  2.88998647 
Li  0.00764420  0.00764342  0.03627239 
H  0.01970090  0.01970037  2.82236966 









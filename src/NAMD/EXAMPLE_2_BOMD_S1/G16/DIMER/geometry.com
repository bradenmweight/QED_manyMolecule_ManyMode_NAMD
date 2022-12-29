%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 1000

0 1
Li  2.05784675  1.82840903  2.28487710 
H  0.27938737  1.85966853  1.71367059 
Li  2.03343844  1.83295548  2.27850755 
H  0.46210587  1.84280791  1.77205154 









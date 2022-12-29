%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/STO-3G nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 1000

0 1
Li  1.75369494  1.69545741  2.35350149 
H  2.37457318  2.77621452  1.24110554 
Li  1.76309086  1.70981190  2.34273130 
H  2.32428623  2.69174964  1.32981945 









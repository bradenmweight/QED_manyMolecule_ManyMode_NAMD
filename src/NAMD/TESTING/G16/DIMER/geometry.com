%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 11

0 1
Li  0.00919629  0.00919615  0.05058556 
H  0.02347538  0.02347770  2.73818950 
Li  0.01013303  0.01013287  0.06032854 
H  0.02570501  0.02570754  2.67971987 









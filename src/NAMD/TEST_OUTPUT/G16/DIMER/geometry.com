%chk=geometry.chk
%nprocshared=1
%mem=5GB

# B3LYP/6-31G* nosymm iop(2/12=3,3/33=1) guess=only pop=full

MD Step 123

0 1
Li  1.02091832  1.02057256  0.81850136 
F  0.61266711  0.61275768  4.68667607 
Li  1.02988743  1.02953820  0.81916510 
F  0.61746932  0.61756086  4.69451332 









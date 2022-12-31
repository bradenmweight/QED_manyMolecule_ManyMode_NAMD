%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# B3LYP/6-31G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 122

0 1
Li  1.02091832  1.02057256  0.81850136 
F  0.61266711  0.61275768  4.68667607 









%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.26749684  -0.08072639  -0.03614803 
H  0.05106432  -1.12191401  -0.16431557 
H  -1.37170355  0.15678200  -0.25993580 
C  0.72022969  0.88821521  0.02142758 
H  0.32006224  1.89888203  -0.27237982 
H  1.79210441  0.54610851  -0.21064827 









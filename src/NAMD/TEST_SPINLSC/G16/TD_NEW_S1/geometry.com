%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.27003675  -0.08546663  -0.03439252 
H  0.05807867  -1.11772663  -0.15572286 
H  -1.36365347  0.17588947  -0.27978038 
C  0.72186284  0.89218922  0.02180105 
H  0.32332196  1.90135619  -0.26088867 
H  1.78731277  0.53626633  -0.22657284 








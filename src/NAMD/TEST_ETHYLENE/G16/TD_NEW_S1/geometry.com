%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311+G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  0.98397243  -0.17551555  0.11279196 
H  -8.85548803  -6.38748218  -1.35078773 
H  -7.96296321  -6.15430436  -1.74621927 
C  0.32766879  1.09555253  0.09156225 
H  0.05148747  1.79690361  -0.73501861 
H  1.19832796  -0.17553435  1.46221598 









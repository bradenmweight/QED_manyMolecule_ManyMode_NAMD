%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 999

0 1
C  2.72953842  1.87134338  2.07945838 
H  -1.05560613  -1.02062931  2.50032533 
H  -1.04327178  -0.67903808  1.52510196 
C  1.53942355  3.05223369  2.06007944 
H  1.75350845  3.23019869  3.16005581 
H  1.57091915  3.70494375  1.23458137 









%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311+G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 571

0 1
C  -0.14366031  -0.65356283  0.10338814 
H  0.47763801  -1.60330036  -0.51561816 
H  -1.07006809  -1.30484788  0.29285911 
C  0.08728505  0.86410079  -0.12461114 
H  0.97159988  0.73680191  -0.60247964 
H  0.36486140  -0.28561552  1.14012049 









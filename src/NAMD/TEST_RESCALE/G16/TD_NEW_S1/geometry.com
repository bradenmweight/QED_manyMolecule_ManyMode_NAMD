%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=24
%mem=5GB

# WB97XD/6-311++G** SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 2

0 1
C  -0.33701235  -0.58693592  0.00072688 
H  0.20204464  -1.50966191  0.00620281 
H  -1.40221711  -0.58137728  0.00621080 
C  0.33847282  0.58839289  0.00072523 
H  -0.18962328  1.52203606  0.00620827 
H  1.41451752  0.59377405  0.00621607 









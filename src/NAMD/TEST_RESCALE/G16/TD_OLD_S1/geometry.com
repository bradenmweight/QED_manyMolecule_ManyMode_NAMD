%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=24
%mem=5GB

# WB97XD/6-311++G** SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1

0 1
C  -0.33654406  -0.58641310  0.00112017 
H  0.19946486  -1.51149325  0.00386092 
H  -1.40395661  -0.58364762  0.00386485 
C  0.33878283  0.58865206  0.00111934 
H  -0.19174038  1.51922151  0.00386372 
H  1.41168145  0.59137387  0.00386764 








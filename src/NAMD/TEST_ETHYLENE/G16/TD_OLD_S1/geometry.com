%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311+G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 570

0 1
C  -0.14135835  -0.65221934  0.10343543 
H  0.50404545  -1.59902827  -0.52219943 
H  -1.09845734  -1.28571447  0.29485283 
C  0.08333815  0.86612546  -0.12572233 
H  0.97684639  0.69912124  -0.59844541 
H  0.38120638  -0.31141348  1.15337046 









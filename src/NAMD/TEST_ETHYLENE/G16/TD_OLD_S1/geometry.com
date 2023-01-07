%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# WB97XD/6-311+G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 999

0 1
C  0.97779819  -0.17726116  0.11154703 
H  -8.82025202  -6.35565316  -1.33796445 
H  -7.94757265  -6.15416947  -1.75374135 
C  0.33232557  1.09186052  0.09587197 
H  0.02453481  1.81000102  -0.76594260 
H  1.19275317  -0.15579542  1.45134437 









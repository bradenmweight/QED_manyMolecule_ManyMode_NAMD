%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# B3LYP/6-31G* SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 123

0 1
Li  1.02988743  1.02953820  0.81916510 
F  0.61746932  0.61756086  4.69451332 









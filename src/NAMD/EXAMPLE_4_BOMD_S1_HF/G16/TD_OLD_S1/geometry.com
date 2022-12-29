%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/3-21G* SCF=XQC TD=(singlets,nstates=8,root=1) FORCE nosym pop=full guess=read

MD Step 999

0 1
Li  2.21337501  2.21008278  1.69066196 
F  0.66499332  0.66612330  3.85590398 









%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 11

0 1
Li  0.01013303  0.01013287  0.06032854 
H  0.02570501  0.02570754  2.67971987 









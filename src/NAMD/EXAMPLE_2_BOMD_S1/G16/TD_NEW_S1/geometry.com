%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1000

0 1
Li  2.03343844  1.83295548  2.27850755 
H  0.46210587  1.84280791  1.77205154 









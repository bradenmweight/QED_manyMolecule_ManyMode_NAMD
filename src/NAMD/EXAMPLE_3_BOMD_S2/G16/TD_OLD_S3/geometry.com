%oldchk=../TD_NEW_S1/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=10GB

# BLYP/STO-3G SCF=XQC TD=(read,singlets,nstates=4,root=3) FORCE nosym pop=full guess=read

MD Step 999

0 1
Li  1.75369494  1.69545741  2.35350149 
H  2.37457318  2.77621452  1.24110554 








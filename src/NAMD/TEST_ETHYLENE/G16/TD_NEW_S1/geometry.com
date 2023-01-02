%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/3-21G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 158

0 1
C  -0.27400519  -0.62586169  0.14836281 
H  -0.26687808  -1.55201715  -0.59871626 
H  -1.27306886  -0.93896664  0.68494287 
C  0.38475004  0.75933116  -0.22226671 
H  -0.06224675  1.66433480  -0.51163353 
H  0.33040796  -0.70767944  1.36913686 









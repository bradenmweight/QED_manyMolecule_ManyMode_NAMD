%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BLYP/3-21G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 157

0 1
C  -0.27689639  -0.62549364  0.14375261 
H  -0.25883983  -1.55782415  -0.58853549 
H  -1.25224175  -0.93432072  0.70195358 
C  0.38745296  0.75726598  -0.22141223 
H  -0.07819901  1.65681579  -0.48706486 
H  0.31975458  -0.67877529  1.36212535 









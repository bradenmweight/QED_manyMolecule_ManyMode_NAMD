%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=1
%mem=5GB

# BPL/sto-3g TD=(singlets,nstates=4,root=1) FORCE nosym pop=full guess=read

MD Step 1

0 1
C  -1.59610985  -0.22883295  0.00000000 
C  -0.20094985  -0.22883295  0.00000000 
C  0.49658815  0.97891805  0.00000000 
C  -0.20106585  2.18742705  -0.00119900 
C  -1.59589085  2.18734905  -0.00167800 
C  -2.29349185  0.97914305  -0.00068200 
H  -2.14586885  -1.18114995  0.00045000 
H  0.34855815  -1.18134595  0.00131500 
H  1.59626815  0.97899805  0.00063400 
H  0.34913415  3.13957005  -0.00125800 
H  -2.14601285  3.13963005  -0.00263100 
H  -3.39309585  0.97932605  -0.00086200 









%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33584103  -0.58129173  0.00513682 
H  0.21495346  -1.52766573  -0.01515690 
H  -1.41898018  -0.58607689  -0.01522382 
C  0.33290454  0.57633486  0.00404093 
H  -0.16150742  1.56741278  -0.01704501 
H  1.44543039  0.65030592  -0.01702512 








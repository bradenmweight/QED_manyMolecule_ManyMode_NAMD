%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33734775  -0.58377359  0.00743700 
H  0.23206655  -1.51907710  -0.02297254 
H  -1.42335685  -0.57305284  -0.02410498 
C  0.33240411  0.57555335  0.00455949 
H  -0.14837562  1.57427061  -0.02689408 
H  1.44174013  0.65897204  -0.02581042 








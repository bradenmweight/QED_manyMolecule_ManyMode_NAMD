%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33592478  -0.58103950  0.00512533 
H  0.21529890  -1.52980188  -0.01553388 
H  -1.42110725  -0.59327329  -0.01622043 
C  0.33335208  0.57811439  0.00416565 
H  -0.16219651  1.56217564  -0.01774967 
H  1.44274822  0.63984288  -0.01712090 








%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33631953  -0.58209307  0.00578475 
H  0.21959139  -1.52829533  -0.01704588 
H  -1.42199760  -0.58587493  -0.01777792 
C  0.33374069  0.57795447  0.00455838 
H  -0.16105064  1.56727409  -0.01980102 
H  1.44388648  0.64590729  -0.01893332 








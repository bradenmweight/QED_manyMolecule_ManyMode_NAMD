%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33792504  -0.58385970  0.00576969 
H  0.23207414  -1.52621926  -0.01933067 
H  -1.42543088  -0.57708805  -0.02008336 
C  0.33255658  0.57705371  0.00411712 
H  -0.14779135  1.56350845  -0.02323417 
H  1.43802621  0.65379812  -0.02226506 








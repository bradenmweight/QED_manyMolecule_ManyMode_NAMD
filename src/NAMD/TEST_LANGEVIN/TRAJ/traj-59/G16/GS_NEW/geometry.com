%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33496706  -0.57994625  0.00378721 
H  0.20373534  -1.53512657  -0.00942023 
H  -1.42177295  -0.60121596  -0.00993633 
C  0.33407587  0.57878220  0.00304626 
H  -0.17591553  1.55995193  -0.01050509 
H  1.44606445  0.63175452  -0.01006283 








%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33533907  -0.58025548  0.00728990 
H  0.21414942  -1.53515806  -0.02582179 
H  -1.42645426  -0.59850375  -0.02815072 
C  0.33393090  0.57915964  0.00614634 
H  -0.16665107  1.56071937  -0.03077218 
H  1.44270305  0.63298021  -0.02838384 








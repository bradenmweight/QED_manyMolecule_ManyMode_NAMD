%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33525985  -0.58009644  0.00609371 
H  0.21350467  -1.53234434  -0.01732457 
H  -1.42243104  -0.59098571  -0.01985216 
C  0.33248038  0.57638603  0.00407086 
H  -0.16307548  1.56186781  -0.02346263 
H  1.44445865  0.64501964  -0.02113173 








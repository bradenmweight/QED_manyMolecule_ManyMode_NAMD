%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33526904  -0.57974008  0.00354707 
H  0.20996495  -1.53332530  -0.00796287 
H  -1.41867249  -0.59553265  -0.00843299 
C  0.33258897  0.57682494  0.00281809 
H  -0.16578798  1.56336484  -0.00935825 
H  1.44760175  0.64139588  -0.00892628 








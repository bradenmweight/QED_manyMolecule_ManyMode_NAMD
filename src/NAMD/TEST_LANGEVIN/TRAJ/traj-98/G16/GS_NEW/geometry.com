%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33527877  -0.58046729  0.00973173 
H  0.20806193  -1.53398077  -0.03964322 
H  -1.42113989  -0.59764893  -0.04141112 
C  0.33425460  0.57931065  0.00837092 
H  -0.16952594  1.55860520  -0.04447130 
H  1.44269193  0.63469466  -0.04229236 








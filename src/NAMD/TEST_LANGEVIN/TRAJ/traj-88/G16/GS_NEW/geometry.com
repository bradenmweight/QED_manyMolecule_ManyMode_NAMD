%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33635446  -0.58103805  0.00784256 
H  0.21932989  -1.53571488  -0.03054373 
H  -1.42316341  -0.58911216  -0.03075806 
C  0.33335275  0.57817422  0.00605383 
H  -0.15838719  1.55552323  -0.03251087 
H  1.43769846  0.64313880  -0.03206644 








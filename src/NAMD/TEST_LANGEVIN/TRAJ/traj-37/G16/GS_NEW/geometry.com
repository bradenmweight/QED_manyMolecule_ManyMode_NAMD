%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33693955  -0.58316810  0.00267187 
H  0.22142100  -1.52273760  -0.00161443 
H  -1.41923480  -0.58210125  -0.00215483 
C  0.33284857  0.57707773  0.00199374 
H  -0.15855021  1.57110101  -0.00435084 
H  1.44898569  0.65017973  -0.00360200 








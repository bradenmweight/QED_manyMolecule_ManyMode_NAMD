%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33694705  -0.58274307  0.00698462 
H  0.22654983  -1.52259999  -0.02240655 
H  -1.42080068  -0.58132541  -0.02398165 
C  0.33260916  0.57692322  0.00517103 
H  -0.15054270  1.56753878  -0.02674673 
H  1.44303453  0.65228604  -0.02513509 








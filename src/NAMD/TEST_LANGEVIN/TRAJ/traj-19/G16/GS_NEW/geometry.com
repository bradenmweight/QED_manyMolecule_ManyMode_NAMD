%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33701524  -0.58267008  0.00716946 
H  0.22658574  -1.52369876  -0.02395004 
H  -1.42223036  -0.58462830  -0.02611863 
C  0.33298800  0.57764983  0.00533714 
H  -0.15214073  1.56619973  -0.02833362 
H  1.44025835  0.64643605  -0.02614565 








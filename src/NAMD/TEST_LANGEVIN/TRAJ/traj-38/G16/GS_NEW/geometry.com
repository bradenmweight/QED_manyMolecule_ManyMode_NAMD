%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33516947  -0.58026383  0.00595075 
H  0.20576048  -1.53319661  -0.02081317 
H  -1.42199716  -0.60059913  -0.02195370 
C  0.33420884  0.57901670  0.00483500 
H  -0.17201584  1.55835090  -0.02262804 
H  1.44150844  0.63210825  -0.02131325 








%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33572637  -0.58034596  0.00553091 
H  0.21592770  -1.53353278  -0.01783038 
H  -1.42368573  -0.59506815  -0.01886982 
C  0.33319805  0.57764492  0.00463305 
H  -0.16278645  1.56391312  -0.02114218 
H  1.44391577  0.64010892  -0.02004631 








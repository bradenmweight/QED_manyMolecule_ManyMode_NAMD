%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=2
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33705979  -0.58307979  0.00421428 
H  0.22410823  -1.51988114  -0.00920705 
H  -1.41863541  -0.58058769  -0.01010689 
C  0.33223195  0.57571529  0.00284417 
H  -0.15336921  1.57502307  -0.01219633 
H  1.44670237  0.65448316  -0.01131462 








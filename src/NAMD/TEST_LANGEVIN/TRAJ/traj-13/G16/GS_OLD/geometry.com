%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33623161  -0.58168141  0.01004920 
H  0.21746505  -1.52221896  -0.03589395 
H  -1.41547722  -0.58522750  -0.03805787 
C  0.33229219  0.57592424  0.00658013 
H  -0.15846539  1.57002752  -0.04139378 
H  1.44724276  0.64985758  -0.03896409 








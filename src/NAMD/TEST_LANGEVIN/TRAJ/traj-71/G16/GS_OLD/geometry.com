%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33510568  -0.57972221  0.00687156 
H  0.21410175  -1.53363639  -0.02316346 
H  -1.42471060  -0.59850975  -0.02545086 
C  0.33322102  0.57794363  0.00572432 
H  -0.16531847  1.56498021  -0.02827798 
H  1.44561832  0.63560517  -0.02595667 








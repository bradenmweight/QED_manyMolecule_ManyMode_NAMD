%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33579806  -0.58113315  0.00816161 
H  0.21253131  -1.53137446  -0.03162177 
H  -1.42172320  -0.59327424  -0.03282590 
C  0.33350367  0.57813529  0.00646907 
H  -0.16413224  1.56134297  -0.03490459 
H  1.44227609  0.64064951  -0.03336692 








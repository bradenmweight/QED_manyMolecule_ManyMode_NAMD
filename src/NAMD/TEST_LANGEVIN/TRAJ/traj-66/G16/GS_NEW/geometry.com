%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33416032  -0.57912761  0.00816846 
H  0.19125071  -1.54002683  -0.03464579 
H  -1.41778321  -0.60282013  -0.03505750 
C  0.33435032  0.57850550  0.00678792 
H  -0.18268061  1.55348031  -0.03618769 
H  1.44365411  0.63348451  -0.03562717 








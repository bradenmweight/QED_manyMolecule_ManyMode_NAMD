%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33575280  -0.58083586  0.00917851 
H  0.22145312  -1.52828805  -0.03300337 
H  -1.42405224  -0.58873454  -0.03511951 
C  0.33273694  0.57677363  0.00659474 
H  -0.15758306  1.56677290  -0.03826011 
H  1.44176066  0.64429226  -0.03592131 








%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33877210  -0.58638413  0.00611988 
H  0.23523102  -1.51063686  -0.01693846 
H  -1.41729680  -0.56447801  -0.01846631 
C  0.33237390  0.57603307  0.00370321 
H  -0.14231186  1.57333260  -0.02150477 
H  1.44087218  0.66536867  -0.01988318 








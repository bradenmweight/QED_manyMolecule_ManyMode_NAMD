%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33723641  -0.58345877  0.00735012 
H  0.22846580  -1.51728788  -0.02384377 
H  -1.41927555  -0.57568346  -0.02500429 
C  0.33178716  0.57509235  0.00462410 
H  -0.14875049  1.57429685  -0.02766640 
H  1.44427087  0.65814301  -0.02639602 








%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  -0.33643962  -0.58253097  0.00656171 
H  0.22025206  -1.52797786  -0.02200365 
H  -1.42297753  -0.58102724  -0.02266816 
C  0.33334400  0.57723290  0.00482187 
H  -0.15840512  1.56269622  -0.02536372 
H  1.43891589  0.65034252  -0.02470976 








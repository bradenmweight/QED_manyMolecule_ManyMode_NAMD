%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33639553  -0.58238826  0.00597403 
H  0.21909738  -1.52549044  -0.01877946 
H  -1.42034673  -0.58233550  -0.01912900 
C  0.33280723  0.57635254  0.00426817 
H  -0.15899077  1.57036049  -0.02067065 
H  1.44617929  0.65257325  -0.02027545 








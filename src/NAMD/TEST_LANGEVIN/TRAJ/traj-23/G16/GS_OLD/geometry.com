%oldchk=../GS_OLD/geometry.chk
%chk=geometry.chk
%nprocshared=4
%mem=5GB

# WB97XD/6-311G SCF=XQC FORCE nosym pop=full guess=read

MD Step 999

0 1
C  -0.33685233  -0.58153398  0.00690250 
H  0.22205985  -1.52910069  -0.02625244 
H  -1.42187339  -0.59033527  -0.02822576 
C  0.33183661  0.57635038  0.00520115 
H  -0.15462495  1.56553193  -0.03100703 
H  1.44424037  0.64570728  -0.02870737 








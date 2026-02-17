#############################################################################
### glaxo_bases.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
define $B1 [NH2][CX4]
define $B2 [NH]([CX4])[CX4]
define $B3 [NX3]([CX4])([CX4])[CX4]
define $B4SUB [C,c](=N)N
define $B5EXC [C,c](=N)N[C,S]=O
define $B6 [$B4SUB;!$B5EXC]
define $B7EXC n(:c)(:c):a
define $B7 [nH0;!$(n-C);!$B7EXC]1ccccc1
define $B8EXC [N,n;+1]
define $B9 [$([NH2]!:c),$([NH1]([CX4])!:c),$([NH0]([CX4])([CX4])!:c)]
[$B1,$B2,$B3,$B6,$B7,$B9;!B8EXC]

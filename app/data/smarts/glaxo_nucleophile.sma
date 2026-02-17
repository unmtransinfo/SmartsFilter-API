#############################################################################
### glaxo_nucleophile.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
define $N1 [NH2][CX4]
define $N2 [NH]([CX4])[CX4]
define $N4EXC N-[O,C,N,S]
define $N5EXC N-[C,c,N]=[C,c,N,n,O,S]
define $N6 [OH1][C,c,N;!$(C=O)]
define $N7EXC [OH1]C=C
define $N8EXC [OH1]NC=[O,S]
define $N9 [$([NH2]!:c),$([NH1]([CX4])!:c),$([NH0]([CX4])([CX4])!:c)]

# alcohols, hydroxylamines but excludes e.g., carboxylic acids
[$N1,$N2,$N6,$N9;!$N4EXC;!$N5EXC;!$N7EXC;!$N8EXC]

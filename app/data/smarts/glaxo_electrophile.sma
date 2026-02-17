#############################################################################
### glaxo_electrophile.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
# alkyl and aryl ketones and aldehydes
define $E1 [C;H1](=[O,S])[C,c]
# e.g, carboxylic esters
define $E2 [C,P;H1](=[O,S])[O,S]
# e.g., carbonates
define $E3 C(=O)([C,c,O,S])[C,c,O,S]
define $E4EXC C(=O)[OH1]
define $E5EXC C(=O)[SH1]
define $E6 C(=[O,S])(N)Oc
# e.g., aryl carbamates
define $E7 C1(=O)NS(=O)(=O)[C,c]=,:[C,c]1
define $E8SUB P(=O)[O,S]
define $E9EXC P[OH1]
define $E10 [$E8SUB;!$E9EXC]
define $E11 c(=O)(~c)~c
define $E12EXC [$(c1(=O)ccn([C,c])cc1),$(c1(=O)n([C,c])cccc1)]
define $E12 [$E11;!$E12EXC]
# imides
define $E13 C(=O)-N-C=O
[$E1,$E2,$E3,$E6,$E7,$E10,$E12,$E13;!$E4EXC;!$E5EXC]

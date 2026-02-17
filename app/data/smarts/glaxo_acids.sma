#############################################################################
### glaxo_acids.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
[OH1][P,C,S]=O A1, e.g. carboxylic acid
[NH1]([P,S]=O)[P,S]=O A2, e.g. imide
[nH]1cnoc1=O A3
[OH1]C1=NC=NO1 A4
[NH1]1C=NOS1=O A5
[OH1]C1=NC(=O)CC1=O A6
[OH1]C1NC(=O)C(=O)C1 A7
[nH1]1ncoc1=O A8
[OH1]C1=NN=CO1 A9
[nH1]1[nH]cnc1=O A10
[OH1]C1=N[NH1]C=N1 A11
[OH1]C1=NOC=C1 A12
[nH1]1occc1=O A13
[OH1]c1oncc1 A14
[nH1]1ccc(=O)o1 A15
[nH1]nnn A16 tetrazole
[nH1](n)nn A17
[OH1]C1=NC(=O)NO1 A18
[OH1]C1=NC(=O)ON1 A19
[nH1]1cnnc1C(F)(F)F A20
[nH1]1cnc(n1)C(F)(F)F A21
[nH1]1C(=O)CC(=O)O1 A22
[OH1]C1=CC(=O)NO1 A23
[OH1]C1=CC(=O)ON1 A24
[NH1]1C(=O)c2ccccc2S1(=O)=O A25 benzosulphimide
[OH1]C1=NS(=O)(=O)c2ccccc21 A26
[OH1]C1=NC(=O)c2ccccc21 A27
[OH1]C1=COC=CC1=O A28
[OH1]C1=NSN=C1 A29
[OH1]NC=O A30 hydroxamic acid
[NH]S(=O)(=O)C(F)(F)F A31 trifluoromethyl sulphonamide
[NH](c)S(=O)=O A32 aryl sulphonamide
[OH1]c1c[c,n]ccc1 A33 phenol

#############################################################################
### glaxo_reactive.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
[Br,Cl,I][CX4;H1,H2] R1 reactive alkyl halides
[S,C](=[O,S])[F,Br,Cl,I] R2 acid halides
O=CN=[N+]=[N-] R3 carbazides
COS(=O)O[#6] R4 sulphate esters
COS(=O)(=O)[#6] R5 sulphonates
C(=O)OC=O R6 acid anhydrides
OO R7 peroxides
C(=O)Oc1c(F)c(F)c(F)c(F)c1F R8 pentafluorophenyl esters
C(=O)Oc1ccc(N(=O)=O)cc1 R9 paranitrophenyl esters
C(=O)Onnn R10 esters of HOBT
N=C=[S,O] R11 isocyanates & isothiocyanates
OS(=O)(=O)C(F)(F)F R12 triflates
P(=S)(S)S R13 Lawesson's reagent and derivatives
NP(=O)(N)N R14 phosphoramides
cN=[N+]=[N-] R15 aromatic azides
C(=O)C[#7+] R16 beta carbonyl quaternary nitrogen
[N&R0][N&R0]C=O R17 acylhydrazide
[C+,Cl+,I+,P+,S+] R18 quaternary C,Cl,I, P or S
C=P R19 phosphoranes
[Cl]C([C&R0])=N R20 chloramidines
[ND2]=O R21 nitroso
[P,S][Cl,Br,F,I] R22 P/S halides
N=C=N R23 carbodimide
[N+]#[C-] R24 isonitrile
C(=O)N(C(=O))OC=O R25 triacyloximes
N#CC[OH] R26 cyanohydrins
N#CC=O R27 acyl cyanides
S(=O)(=O)C#N R28 sulfonyl cyanides
P(OCC)(OCC)(=O)C#N R29 cyanophosphonates
[N&R0]=[N&R0]C#N R30 azocyanamides
[N&R0]=[N&R0]CC=O R31 azoalkanals

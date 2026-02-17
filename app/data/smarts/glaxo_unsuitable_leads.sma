#############################################################################
### glaxo_unsuitable_leads.sma - 
###
### ref: "Strategic Pooling of Compounds for High-Throughput Screening",
### Mike Hahn et al., J. Chem. Inf. Comput. Sci. 1999, 39, 897-902
#############################################################################
[CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0] I1 aliphatic methylene chains 7 or more long
[C,S,P](=O)[OH].[C,S,P](=O)[OH].[C,S,P](=O)[OH].[C,S,P](=O)[OH] I2 compounds with 4 or more acidic groups
[O;R1][C;R1][C;R1][O;R1][C;R1][C;R1][O;R1] I3 crown ethers
SS I4 disulphides
[SH] I5 thiols
C1[O,S,N]C1 I6 epoxides, thioepoxides, aziridines
c([OH])c([OH])c[OH] I7 2,4,5 trihydroxyphenyl
c([OH])c([OH])cc[OH] I8 2,3,4 trihydroxyphenyl
N=NC(=S)N I9 hydrazothiourea
SC#N I10 thiocyanate
cC[N+] I11 benzylic quaternary nitrogen
C[O,S;R0][C;R0]=S I12 thioesters
N[CH2]C#N I13 cyanamides
C1(=O)OCC1 I14 four membered lactones
P(=O)([OH])OP(=O)[OH] I15 di and triphosphates
N1CCC1=O I16 betalactams

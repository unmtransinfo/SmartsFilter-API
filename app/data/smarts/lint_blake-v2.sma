#############################################################################
### lint_blake.sma - smarts translated from "lint_sln.spl" by James Blake.
### These smarts match compounds considered undesirable in a HTS.
###
### Jeremy Yang
###   1 Feb 2008
#############################################################################
O~N(=O)-c(:*):* aromatic NO2
[2H] deuterium
[13#6] C13
n1c(cccc1)Cl 2-chloropyridine
[NH2D1]-c(:*):* aniline
[Si,B,Se] Si,B,Se atoms
[!#6]-[CH2]-N1C(=O)CCC(=O)1 hetero imides
O[CH2][CH2]O-!@[CH2][CH2]O poly ethers
C-!@[NX2]=[C!R,#1]-C acyclic imines
[S,P](=O)OC alkyl esters of S or P
P(=[O,S])[C,N]([C,N])[C,N] ugly P compounds
[N;!$(N-[C,S]=*)]-,=;!@[N;!$(N-[C,S]=*)] acyclic N-,=N and not N bound to carbonyl or sulfone
N-!@[CX4]-!@N acyclic N-C-N
N-!@[SX2]-* acyclic N-S
[N,S,O][CH2][CH2]-[F,Cl,Br,I] mustards
O=[C!H0] aldehyde
O=[CX3]-!@[CX3]=O 1,2-dicarbonyl not in ring
NC(OC([CH3])([CH3])[CH3])=O carbamate, T-boc Protected
NC(O[CH2]c1ccccc1)=O carbamate, CBZ Protected
OC(=O)-!@[NX3] carbamate include di-substitued N
O=C-[F,Cl,Br,I] acyl halide
[CH2]-[Cl,Br,I] alkyl halide
O=C-C-[F,Cl,Br,I] alpha halo carbonyl
S(=O)(=O)-[F,Cl,Br,I] sufonyl halide
[ND1]=C-!@[SX2]-[CH2D2] N:C-SCH2 groups
N-C(=S)-N
[CH2D1]=[CD2]-!@* terminal vinyl
S=P~*
S-C#N thio cyanates
*-[S!H0] thiols
[Sv4](=O)(-!@[!#1])-!@[!#1] thionyl
N-[F,Cl,Br,I] n-haloamines
N-C-[F,Cl,Br,I,$(C#N)] N-C-Hal or cyano methyl
[$(C#N),$(N(~O)~O),$(C=O),$(S(=O)=O),$(C(F)(F)F),Cl][C!H0]=[C!H0] alpha beta-unsaturated ketones; center of Michael reactivity
[C;!$(C=*)][C!R](=O)[CH2D2] aliphatic ketone not ring and not di-carbonyl
C[C!R](=O)[O!R][CH2D2] aliphatic ester, not lactones
[CH2][CH2][CH2]-!@[CH2][CH2][CH2] long aliphatic chain, 6+
O=[#6]1[#6]:,=[#6][#6](=O)[#6]:,=[#6]1 quinones
C=C-!@O-* acyclic C=C-O
[NH1X2]=[C!R;!$(C(-N)(=[NH1])-N)] acyclic C=N-H
O-!@[N;!$(N(=O)=O);!$([N+](=O)[O-])] acyclic NO not nitro
S-S
O~O
C-C(=O)[SD2] thioester
N~1~*~*1 aziridine-like N in 3-membered ring
O~1C~*1 epoxides
I-c:* aryl iodide
Br-c:* aryl bromide
a1aaa2a(a1)aaa(a2):a multiple aromatic rings
[!#1]-,:1-:a2a(-:[!#1]:3:a1aaaa3)aaaa2 multiple aromatic rings
[P,S](~O)(~O)~O S/PO3 groups
C12CC3CC(C1)CC(C2)C3 adamantyl
C#N.C#N too many cyano Groups (>1)
[CX3](=O)[OH1].[CX3](=O)[OH1] too many COOH groups (>1)
[NH2][CX4]C(=O)O amino acid
Cl~O chlorates
[F,Cl,Br,I].[F,Cl,Br,I].[F,Cl,Br,I].[F,Cl,Br,I] high halogen content (>3)

#############################################################################
### oprea_filters.sma
###
### See:
### http://hatch.health.unm.edu/dbcwiki/index.php/Roadrunner
###
#############################################################################
#
# Nitro both forms
#
define $no2A [N;+1](=O)[O;-1]
define $no2B N(=O)=O
define $no2 [$no2A,$no2B]
#
# amines
#
define $spone [S,P]=O
define $cene C=[O,S,C]
define $elneg [$cene,$spone]
define $amine1 [N;D1;H3&+1,H2&+0][C,c;!$elneg]
define $amine2 [N;D2;H2&+1,H&+0]([C,c;!$elneg])[C,c;!$elneg]
define $amine3 [N;D3;H&+1,H0&+0]([C,c;!$elneg])([C,c;!$elneg])[C,c;!$elneg]
#
# not carbon
#
define $notc [!#6]
#
#   Al/Ar Amine 
#
define $cnh [$amine1,$amine2,$amine3]
define $cxnh [$notc;!$cnh]~[C,c]
#
#   Halogen 
#
define $nothal [!Cl;!Br;!I;!F]
define $anyhal [Cl,Br,I,F]
define $halnotF [Cl,Br,I]
#
#   Cyanide 
#
define $cyanide C#N
#
#   CH2Hal and CH(Hal)2 and C(Hal)3 
#
define $ewx1 [C;H2][$anyhal]
define $ewx2 [C;H1]([$anyhal])[$anyhal]
define $ewx3 C([$anyhal])([$anyhal])[$anyhal]
define $ewcl1 [C;X4][$halnotF]
define $ewcl2 [C;X4]([$halnotF])[$halnotF]
define $ewcl3 [C;X4]([$halnotF])([$halnotF])[$halnotF]
#
#   SO3C(Hal)3 SO3CH(Hal)2 SO3CH2Hal SO3CH3 
#
define $smarty1 S(=O)(=O)(O)C([$anyhal])([$anyhal])[$anyhal]
define $smarty2 S(=O)(=O)(O)[C;H1]([$anyhal])[$anyhal]
define $smarty3 S(=O)(=O)(O)[C;H2][$anyhal]
define $smarty4 S(=O)(=O)(O)[C;H3]
#
#
(*).(*) disconnected structures
#
#   includes salts, hydrates, racemates 
#
#   isotopes 
#
define $iso1 [2H]
define $iso2 [3H]
define $iso3 [13C]
define $iso4 [13c]
define $iso5 [14C]
define $iso6 [14c]
define $iso7 [15N]
define $iso8 [15n]
[$iso1,$iso2,$iso3,$iso4,$iso5,$iso6,$iso7,$iso8] isotopes
#
#   unacceptable elements
#
define $AtmOK [C,c,N,n,O,o,F,Si,P,S,s,Cl,Br,I,#1]
#
[A,a;!$AtmOK] AtmCrp
#
#   heteroatoms = polar_atoms N, O, P, S
#
define $polatom [#7,#8,#15,#16]
#
### THIS CANNOT BE CORRECT!
#	[#6]~[$polatom] cpol
#
#   quaternary atoms; excludes sulfoxides 
#
[+1;#6,F,Cl,Br,I,#8,#16;!$([S+][O-])] poschatom 
#
#   compounds with >=4 acidic groups 
#
[C,S,P](=O)[OH].[C,S,P](=O)[OH].[C,S,P](=O)[OH].[C,S,P](=O)[OH] more4acid
#
#
#   hydroxy OH/SH/ethers (bin if <=1 and no other het ie if cor<=1 &&
#    cxor=0)
#   what about if any halogens in addition to 1*O/S? 
#
define $acid1 [S,O;H1]C(=[S,O])[C,c] 
define $acid2 [S,O;H0;-1]C(=[S,O])[C,c] 
define $ester [S,O]([C,c])C(=[S,O])[C,c,N] 
define $coh [O,S;X2;H;!$acid1][C,c] 
define $cxoh [$notc;!$coh]~[C,c] 
define $cor [O,o,S,s;X2;!$ester]([C,c])[C,c] 
define $cxor [$notc;!$cor]~[C,c]
#
#   2. quaternary (bin if <=1 and no other het ie if cquatern<=1 &&
#      cxquater=0)
#   3. or if the only other het is a halogen ie if cquatern<=1 && cxquathl=0 
#
define $quat1 [+1;N;H0;D3]([C,c])[C,c]
define $quat2 [+1;N;D4]([C,c])[C,c]
define $quat3 [+1;n;D3](c)c
define $quat4 [+1;n;D4](c)c
define $noxide1 [n;X3;+1]([O;X1;-1])(c)c
define $noxide2 [N;X4,X3;+1]([O;X1;-1])([C,c])[C,c]
define $cquatern [$quat1,$quat2,$quat3,$quat4;!$noxide1;!$noxide2]
define $cxquater [$notc;!$cquatern;!$anyhal]~[C,c]
#
#   only hetero is 1 COOH or 1 COOR (this def gives count 2 = 1
#      acid/ester)
#   bin if cacid<=2 and cxacid=0 
#
define $acid1a [S,O]=C([S,O;H1])[C,c]
define $acid2a [S,O]=C([S,O;H0;-1])[C,c]
define $estera [S,O]=C([S,O][C,c])[C,c]
define $acidoo [$acid1,$acid2,$acid1a,$acid2a]
define $esteroo [$ester,$estera]
define $cacid [$acidoo,$esteroo]
define $cxacid [$notc;!$cacid]
#
#
#   Garbage, Nitrogen
#
#   hydrazine RRN-NRR not in ring, nor diketo, nor keto/edring, nor
#      diedring 
#
define $hydraz1 [N;X3;R0;!$no2]-[N;X3;!$no2]
define $edring1 cn
define $edring2 ccn
define $edring3 cccn
define $carbonyl C=[O,S]
define $edring [$edring1,$edring2,$edring3]
#
define $hyzinok [N;X3;R0]([$edring,$carbonyl])-[N;X3][$edring,$carbonyl]
define $hydrazin [$hydraz1;!$hyzinok]
[$hydrazin] hydrazin
#
#   azo N=N or diazoniu N#N
#   includes N=NN, N=[N;+1]=[N;-1], [N;-1][N;+1]#N, N=[N:+1][O:-1] 
#
define $azodiaz [N;R0]=,#N
[$azodiaz] azodiaz
#
#   carbodiimide 
#
define $carbdiim N=[C;R0]=N
[$carbdiim] carbdiim
#
#   biguanide 
#
define $biguan N(C(=N)N)C(=N)N
[$biguan] biguan
#
#   isocyanide = isonitril -NC 
#
define $isocn1 N=[C;H0;D1]
define $isocn2 [N;+;X2]#[C;-;H0;D1]
define $isocn [$isocn1,$isocn2]
[$isocn] isocyanide
#
# cyanohydrins and acylcyanides 

define $cyanohydrins N#CC[OH]
define $acylcyanides N#CC=O
#
[$cyanohydrins,$acylcyanides] cyanohydrins and acylcyanides 
#
# nitroso 
#
[N&D2]=O nitroso
#
#  cloramidine 
#
[Cl]C([C&R0])=N cloramidine 
#
#  Schiffs base/imine RN=CRR remove them bin if >=1 (hydrolyse to
#     ketone)
#  (excluding benzodiazepine - probably too narrow an exclusion)
#  excluding cyclic imines
#  excluding ArN=C(Ar) NO - Brian M & Richard S
#
define $bzdiazep N1=Cc2ccccc2NC(=O)C1
define $imine1 [N;X2;H1]=[C;X3]([C,c])[C,c]
define $imine2 [N;X2;R0]([C,c])=[C;H;X3]([C,c])
define $imine3 [N;X2]([C,c])=[C;H2;X3]
define $imine4 [N;X2;R0]([C,c])=[C;X3]([C,c])[C,c]
define $imine5 [N;X2](c)=[C;X3](c)[C,c]
define $imine6 [N;X2](c)=[C;H;X3](c)
define $imine [$imine1,$imine2,$imine3,$imine4;!$imine5;!$imine6]
#
[$imine] imine
#
#  hydrazone RRC=N-NRR and not in ring (the C can be aromatic)
#  hydrazone [N;X2;R0](=[C,c])-[N;X3] (0,$]
#  azine C=N-N=C and azide NNN and not in ring (see azodiaz) 
#
[N;X2;R0](=[C;X3])-[N;X2;R0]=[C;X3] azine
#
#  nitrite 
#
N(=O)-O-[C,c] nitrite
#
#  isourea 
#
[N;X3][C;R0]([O;X2])=[N;X2] isourea
#
#  beta (mesyl or halo) amine or O/S mustard bin if mustard >=1
#  hits on trihal (and dihal) may not be nasty but we still get them 
#
define $mustarda [$halnotF][C;X4][C;X4][$cnh]
define $mustards S(=O)(=O)O[C;X4][C;X4][$cnh]
define $mustard [$mustarda,$mustards]
[$mustard] mustard
#
#  CNOC in chain not in ring - hydroxylamine - remove if NH or OH
#  Hydroxamic acid is OK 
#
define $cnocn1 [N;X3;R0]([C,c])[O;H]
define $cnocn2 [N;X3;R0;H]([C,c])O
define $cnocn3 [N;X3;R0](C=O)O
define $cnocn [$cnocn1,$cnocn2;!$cnocn3]
[$cnocn] cnocn
#
#         Garbage, Halogen
#
#  remove if no halogen >=7
#  more7hal [$anyhal] [7,$]
#  N-halogen & S-halogen & CH2-halogen 
#
define $nhal N[$anyhal]
[$nhal]
define $shal S[$anyhal]
[$shal]
define $phal P[$anyhal]
[$phal]
define $ohal O[$anyhal]
[$ohal]
define $ch2hal [CH2][$anyhal]
[$ch2hal]
#
#  acid halide & thio acid halide 
#
define $thachal [S,O]=C[$anyhal]
[$thachal]
#
#  imine halide 
#
define $iminhal [N;R0]=C[$anyhal]
[$iminhal]
#
#  halopyrimidine ; includes purine bases 
#
define $halopyrimidine [$anyhal]c1ncccn1
[$halopyrimidine]
#
#  sulfonyl halide 
#
define $sulfonhal O=S(=O)[$anyhal]
[$sulfonhal]
#
#  1. triflates SO3CX3 
#
define $trif [$smarty1]
[$trif]
#
#  perhalo ketones 
#
define $haloket [CH2]C(=O)C([$anyhal])([$anyhal])[$anyhal]
[$haloket] perhalo ketones 
#
#  halo methylene ether (nasty) cf Isabel's hlmtheth and hletheth &
#     vinylew
#  halmeo [$ewcl1,$ewcl2,$ewcl3]OC (0,$]
#  alpha halo (hal=1,2,3) ketone bin if >=1 
#
define $halmtha [$ewcl1,$ewcl2,$ewcl3]C=O
define $halmthn [$ewcl1,$ewcl2,$ewcl3]C(=O)N
#
#  halmthk [$halmtha;!$halmthn] (0,$]
#  vinyl and Hal (not nitrile) *** why R0 in Isabel's definitions?
#
define $vinylew [$halnotF][C;R0]=C
[$vinylew]
#
#  CH2Hal: CBr3, CBr2, CBr, CCl2, CCl; hal-C-C-R is ok - OK if >=2*C 
#
define $alkcl3 [C;X4](Cl)(Cl)Cl
define $alkhal1 [C;X4][$halnotF]
define $alkhal2 [C;X4;H2]([$halnotF])[C;X4]
define $alkhal3 [C;X4;H]([$halnotF])([$halnotF])[C;X4]
define $alkhal4 [C;X4]([$halnotF])([$halnotF])([$halnotF])[C;X4]
define $alkhalok [$alkcl3,$alkhal2,$alkhal3,$alkhal4]
#
#  alkhal [$alkhal1;!$alkhalok;!$halmtha;!$halmeo] (0,$]
#  halo (inc F) triazine or 2-, 4-pyrimidine bin all
#  do not bin 2-amino, 4-Cl pyrimidine nor 2-amino, 4-Cl triazine 
#
define $halpyr1 [$anyhal]c1ncncc1
define $halpyn2 [$anyhal]c1nc(N)ncc1
define $halpyn3 [$anyhal]c1ncnc(N)c1
define $halpyn4 [$anyhal]c1ncncc1N
define $halpy1 [$halpyr1;!$halpyn2;!$halpyn3;!$halpyn4]
define $halpyr2 [$anyhal]c1ncccn1
define $halpyn5 [$anyhal]c1nc(N)ccn1
define $halpyn6 [$anyhal]c1ncc(N)cn1
define $halpy2 [$halpyr2;!$halpyn5;!$halpyn6]
define $haltriz1 [$anyhal]c1ncncn1
define $haltrn1 [$anyhal]c1nc(N)ncn1
define $haltri1 [$haltriz1;!$haltrn1]
define $haltriaz [$halpy1,$halpy2,$haltri1]
#
#         Garbage, Oxygen
#
OO peroxides
#
#  Aliphatic C-aldehyde RC(=O)H ; aldehyde is stable if it is cCHO or
#     [N,O,S]CCHO 
#
define $aldnos O=[C;H1]C[N,O,S]
define $aldcx O=[C;H1]C
define $aldehyde [$aldcx;!$aldnos]
[$aldehyde] aldehyde
#
#  anhydride bin if anhydrid>=1 
#
define $anhydride O=C([C,c])OC(=[O,N])[C,c]
[$anhydride] anhydride
#
#  triacyloximes 
#
define $triacyloxime O=CN(C=O)OC(=O)
[$triacyloxime] triacyloxime
#
#  sugar 
#
define $sugar1 O[C;R0][C;R0](O)[C;R0](O)[C;R0](O)[C;R0]-,=O
define $sugar2 O1C([O;R0])C([O;R0])C([O;R0])C([O;R0])C1([O,C;R0])
define $sugar [$sugar1,$sugar2]
[$sugar] sugar
#
#  sugar OCC(O)C(O)C(O)C-,=O (0,$]
#  sugoh [O;X2;H;!$acid1][C,c] *
#  sugor [O;X2;!$ester](C)[C,c] *
#  sugaro [$sugoh,$sugor] (0,$]
#
#  poly (CCO) bin if ccochain>=1 
#
define $cco2 [C;H2][C;H2][O;D2][C;H2][C;H2][O;X2]
define $ccoring [O;D2;R][C;H2][C;H2][O;D2;R][C;H2][C;H2][O;X2;R][C;H2][C;H2]
define $cco5 [O;D2]([$cco2])[$cco2]
define $ccochain [$ccoring,$cco5]
[$ccochain] ccochain
#
#  quinones
#  quinones O=C1[#6]~[#6]C(=O)[#6]~[#6]1 (0,$]
#
#  benzoquinone 
#
define $bquinon4 [N,O;R0]=[C,c]1[C,c]=,:[C,c][C,c](=[C,N,O;R0])[C,c]=,:[C,c]1
define $bquinon2 [N,O;R0]=[C,c]1[C,c](=[C,N,O;R0])[C,c]=,:[C,c][C,c]=,:[C,c]1
define $bquinone [$bquinon4,$bquinon2]
[$bquinone] bquinone
#
#  2,3,4, & 2,4,5 trihydroxyphenyl 
#
define $triOHph1 c([OH])c([OH])c([OH])
define $triOHph2 c([OH])c([OH])cc([OH])
define $triOHph [$triOHph1,$triOHph2]
[$triOHph] triOHph
#
#         Garbage, Sulfur
#
#  SS 
#
define $disulphide [S;$(S[#6,#16]);!$(S(=O))]-[S;$(S[#6,#16]);!$(S(=O))]
[$disulphide] disulphide
#
#  SN 
#
define $SN [S;!$(S(=O));!$(S(=N))]-[N;!$(N=[#6,#7])]
[$SN] SN
#
#  SO 
#
define $SO [S;D1,D2;R0;!$(S(=O))]-[O;D1,D2;R0]
[$SO] SO
#
#  sulfonyl cyanides 
#
define $scn O=S(=O)C#N
[$scn] sulfonyl cyanides 
#
#  sulphonamide imine (Michael) but not sulphonamide amidine nor in ring 
#
define $so2nc1 [S;X4](=O)(=O)N=C
define $so2nc2 [S;X4](=O)(=O)N=CN
define $so2ncr5 [S;X4]1(=O)(=O)N=C[C,c]~[C,c]1
define $so2ncr6 [S;X4]1(=O)(=O)N=C[C,c]~[C,c]~[C,c,N,n]1
define $so2nc [$so2nc1;!$so2nc2;!$so2ncr5;!$so2ncr6]
[$so2nc] so2nc
#
#  isocyanates & isothiocyanates 
#
define $isothcn [S,O]=C=N
[$isothcn] isothcn
#
#  thiocyanate 
#
define $thcn SC#N
[$thcn] thiocyanate 
#
#         Garbage, P(III)
#
#  phosphorous acid/ester 
#
define $pphorous [P;X3]([O;X2])([O;X2])[O;X2]
define $pphonous [P;X3]([O;X2])([O;X2])[C,c]
define $pphinous [P;X3]([O;X2])([C,c])[C,c]
define $pphine [P;X3]([C,c])([C,c])[C,c]
define $phos3 [$pphorous,$pphonous,$pphinous,$pphine]
[$phos3] phosphorous acid/ester 
#
#         Garbage, P(V)
#
#  phosphoric acid/ester 
#
define $pphoric [P;X4](=O)([O;X2])([O;X2])[O;X2]
define $pphonic [P;X4](=O)([O;X2])([O;X2])[C,c]
define $pphinic [P;X4](=O)([O;X2])([C,c])[C,c]
define $pphino [P;X4](=O)([C,c])([C,c])[C,c]
define $pphane [P;X5]([C,c])([C,c])([C,c])([C,c])[C,c]
define $pphonim [P;X4](=O)([O;X2])([N;X3])[C,c]
define $pphnylid [P;X4](=C)([C,c])([C,c])[C,c]
define $phos5 [$pphoric,$pphonic,$pphinic,$pphino,$pphane,$pphonim,$pphnylid]
[$phos5] phosphoric acid/ester 
#
#  Lawesson's reagent and derivatives 
#
S=P(S)S Lawesson's reagent and derivatives 
#
#  phosphoramides 
#
NP(=O)(N)N phosphoramides 
#
#  Disubs phosphine oxide R2P=O 
#
[P;D3](=O)([!O])[!O] Disubs phosphine oxide R2P=O 
#
#  phosphonium cation 
#
[P&+1;D4]([!O])([!O])([!O])[!O] phosphonium cation 
#
#  PS or phosphonium ylide or phosphonium imide 
#
define $pgrot1 P~S
define $pgrot2 P=[N,C]
define $pgrot [$pgrot1,$pgrot2]
[$pgrot] PS or phosphonium ylide or phosphonium imide 
#
#  cyanophosphonates 
#
P(OOC)(OOC)(=O)C#N cyanophosphonates 
#
#         Garbage, Unsaturated
#
#  alpha unsat ketones C=C-C=O Michael acceptors bin if michael>=1
#  if c-C=C-C=O it is possibly stable enough 
#
define $acidlg1 C(=O)[O&H,O&-1]
define $acidlg2 C(=O)[O,S][C,c,N]
define $keto C(=O)([C,c])C
define $ald [C;H](=O)
define $amid C(=O)([N;X3])
define $sone S(=O)(=O)([C,c])[C,O]
define $side S(=O)(=O)([N;X3])
define $michlg1 [$keto,$ald,$amid,$sone,$side,$no2,$acidlg1,$acidlg2]
define $michlg [$keto,$ald,$amid,$sone,$side,$no2,$acidlg1,$acidlg2,$cyanide]
define $michlg2 [$halnotF,$cyanide]
define $michneda [C;H2]
define $michnedb [C;H][!N;!$cor;!$coh;!c]
define $michnedc C([!N;!$cor;!$coh;!c])[!N;!$cor;!$coh;!c]
define $michned [$michneda,$michnedb,$michnedc]
define $michael1 [$michlg]C=[$michned]
[$michael1]
#
define $mich1 [$michlg]C=[$michned]
define $mich2 [$michlg]C=C[$michlg]
define $mich3 [$michlg]C([$michlg])=C
define $michael2 [$mich1,$mich2,$mich3]
[$michael2]
#
#  acetylene 1-,2-cruds 
#
define $michael3 [$michlg]C#C
[$michael3]
#
#  reactive beta keto quaternary O=CCN+ 
#
define $ketoquat [$cquatern]CC=O
[$ketoquat] reactive beta keto quaternary O=CCN+ 
#
#  reactive esters and thioesters RCOOCX2 X=Hal/CN/Ph ?? 
#
define $reactx1 [$halnotF,$cyanide]
define $reactx2 [$anyhal,$cyanide]
define $reactx3 [$anyhal,$cyanide,c]
define $reactst1 O=C[O,S]C[$reactx1]
define $reactst2 O=C[O,S]C([$reactx2])[$reactx2]
define $reactst3 O=C[O,S]C([$reactx3])([$reactx3])[$reactx3]
define $reactst4 O=C(C([$reactx2])[$reactx3])[O,S][C,c]
define $reactst5 O=[C;X3][O,S]c
define $reactst6a O=[C;X3][O,S]C=C
define $reactst6b O=[C;X3]1[O,S]C=CC1
define $reactst6c O=[C;X3]1[O,S]C=CCC1
define $carbamat O=C([O,S])N
define $reactst6 [$reactst6a;!$reactst6b;!$reactst6c]
define $reactest1 [$reactst1,$reactst2,$reactst3,$reactst4,$reactst5,$reactst6]
#
define $reactest [$reactest1;!$carbamat]
[$reactest]
#
#         LARGE RINGS, LARGE ALKYL CHAINS
#
#  Large ring >=C9 (bin if c9str>=1 or c9rng>=1) 
#
define $c9rng1 [r9,r10,r11,r12,r13,r14;C]
define $c9rng [$c9rng1][C;x2][C;x2][C;x2][C;x2][C;x2][C;x2][C;x2][C;x2]
[$c9rng] LARGE RINGS, LARGE ALKYL CHAINS
#
#  crown ethers 
#
define $crownether [O;x2][C;x2][C;x2][O;x2][C;x2][C;x2][O;x2]
[$crownether] crownether
#
#  1. C9 chain Not in any rings = max flex chain 1 descriptor in saSA 
#
define $c9str [C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]
[$c9str] c9str
#
#  straight chain or Al ring no branching (bin if branch=0)
#  straight [A;D2,D1] *
#  branch [!$straight] (0,$]
#  branch [A;D3,D4] (0,$]
#  multi alkene chain CC=CC=CC=C 
#
define $mltalken [C;X3;!R]=C[C;!R]=C[C;!R]=[C;X3]
[$mltalken]
#
#         Ugly Structures
#  annelated rings pyrene, phenanthrene, anthracene, adamantane (bin
#     if >=1) 
#
define $pyrene c1cc2ccc3cccc4ccc(c1)c2c34
define $phenan c12ccccc1ccc3ccccc32
define $anthrac c1c2ccccc2cc3ccccc31
define $big3 c12c3~[N,C,c]~[N,C,c]~[N,C,c]~c1aaaa2aaa3
define $adamant C12CC3CC(C1)CC(C2)C3
define $bigring [$pyrene,$phenan,$anthrac,$big3,$adamant]
[$bigring]
#
#  diazetidinone ring 
#
N1C(=O)NC(=O)1 diazetidinone ring 
#
#  cycloheximide derivatives 
#
O=C1CCCC(N1)=O cycloheximide derivatives 
#
#  epoxides & aziridines & thiiranes & oxazirane 
#
define $epazthi C1[N,O,S]C1
define $oxaz N1[S,O]C1
[$epazthi,$oxaz] epoxides & aziridines & thiiranes & oxazirane 
#
#  4 membered lactones 
#
define $lactone4 O=C1OCC1
[$lactone4] 4 membered lactones 
#
#         Unsuitable Substructures
#
#  diphenyl ethylene cyclohexadiene 
#
define $dyecore1 C1=CCC=CC1=C(c2ccccc2)c3ccccc3
define $dyecore2 C1=CCC=CC1=C(c2ccccc2)C
define $dyecore [$dyecore1,$dyecore2]
[$dyecore]
#
#  p-, p'-dihydroxy biphenyl, stilbene 
#
define $diohbiph [O;H]c1ccc(cc1)c2ccc([O;H])cc2
[$diohbiph]
define $diohstil [O;H]c1ccc(cc1)[C;R0]=[C;R0]c2ccc([O;H])cc2
[$diohstil]
#
#  saponin derivatives 
#
define $saponin O1CCCCC1OC2CCC3CCCCC3C2
[$saponin]
#
#  cytochalasin derivatives 
#
define $cytochalasin O=C1NCC2CCCCC21
[$cytochalasin]
#
#  betalactams 
#
define $blactams O=C1CCN1
[$blactams]
#
#  monesin derivatives 
#
define $monesin O1CCCCC1C2CCCO2
[$monesin]
#
#  squalestatin derivatives 
#
define $squalestatin C12Occc(O1)CC2
[$squalestatin]
#
#  steroid 
#
define $steroid1 [C,c]12~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~1~C(~C3~C~C~2)~C~C~C4CCCC~4~3
define $steroid2 [C,c]12~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~1~C(~C3~C~C~2)~C~C~C4CCCCC~4~3
define $steroid [$steroid1,$steroid2]
[$steroid] steroid
#
#  penicillin or cephalosporin 
#
define $penicil N12C(=O)CC1SC~C2
define $carpenem N12C(=O)CC1CC~C2
define $cephalo N12C(=O)CC1SC~C~C2
define $penceph [$penicil,$carpenem,$cephalo]
[$penceph] penicillin or cephalosporin 
#
#  prostaglandin 
#
define $prostag C1CCC([C;R0]~C~C~C~C~C)C1[C;R0]~C~C~C
[$prostag] prostaglandin 


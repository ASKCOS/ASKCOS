'''
Title:	Extended Functional Groups (EFG): An Efficient Set for Chemical Characterization and Structure-Activity Relationship Studies of Chemical Compounds.
Authors:	Salmina, ES; Haider, N; Tetko, IV;
Journal reference:	Molecules (Basel, Switzerland), 2015; 21 (1);

- Note that NOT and AND must be parsed separately and Hs must be explicit

Manual corrections: 
- "!r0" changed to "!R0"
- 'redux' field added to exclude EFGs that wouldn't impose any cross-reactivity
  e.g., the fact that a molecule is not a hydrocarbon won't make a template not apply
'''

fgroups = [
{	
	'name': 'Hydrocarbons',
	'description': 'Compounds containing only C and H atoms',
	'SMARTS': 'NOT [!#6!#1]',
	'redux': False,
},
{	
	'name': 'Cations',
	'description': 'Any positively charged atom',
	'SMARTS': '[+,++,+3,+4,+5,+6,+7,$([#7v4]),$([#8v3])] AND NOT [$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-]),$([NX2+]#[CX1-]),$([NX2]#[CX1]),$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*+]~[*-])]',
	'redux': False,
},
{	
	'name': 'Anions',
	'description': 'Any negatively charged atom',
	'SMARTS': '[-,--,-3,-4,-5,-6,-7] AND NOT [$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-]),$([NX2+]#[CX1-]),$([NX2]#[CX1]),$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R])),$([*-]~[*+])]',
	'redux': False,
},
{	
	'name': 'Carbonyl compouns: aldehydes or ketones',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[$([H][CX3]([#1,#6])=[OX1]),$([#6&!$(C#N)][CX3](=[OX1])[#6&!$(C#N)])]',
	'redux': True,
},
{	
	'name': 'Aldehydes',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[H][CX3]([#1,#6])=[OX1]',
	'redux': True,
},
{	
	'name': 'Ketones',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$(C#N)][CX3](=[OX1])[#6&!$(C#N)]',
	'redux': True,
},
{	
	'name': 'Thiocarbonyl compounds: thioaldehydes or thioketones',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[$([#1,#6][#6X3](=[SX1])[#1,#6]),$([#6]=[CX2]=[Sv2X1])]',
	'redux': True,
},
{	
	'name': 'Thioaldehydes',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1][CX3]([#1,#6])=[SX1]',
	'redux': True,
},
{	
	'name': 'Thioketones',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$(C#N)][CX3](=[SX1])[#6&!$(C#N)]',
	'redux': True,
},
{	
	'name': 'Imines',
	'description': 'R1, R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[NX2][#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Hydrazones',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[NX2][#7X3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Semicarbazones',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3]([#1,#6])=[NX2][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[CX3](=[OX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Thiosemicarbazones',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl or aromatic carbon',
	'SMARTS': '[#1,#6][CX3]([#1,#6])=[NX2][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[CX3](=[SX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Oximes',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3]([#1,#6])=[NX2][OH1]',
	'redux': True,
},
{	
	'name': 'Oxime ethers',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3]([#1,#6])=[NX2][OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Ketenes',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[OX1]=[CX2]=[CX3]([#1,#6])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Ketene acetal derivatives',
	'description': 'R1, R2 = H, alkyl, aryl; X, Y = any hetero atom',
	'SMARTS': '[#1,#6][CX3]([#1,#6])=[CX3]([!#6!#1])[!#6!#1]',
	'redux': True,
},
{	
	'name': 'Carbonyl hydrates',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([OH1])[OH1]',
	'redux': True,
},
{	
	'name': 'Hemiacetals / Hemiketals',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([OH1])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Acetals / Ketals',
	'description': 'R1, R2 = H, alkyl, aryl; R3, R4 = alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([OX2][#6&!$([CX3]=[OX1,SX1,NX2])])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Hemiaminals',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([OX2][#1,#6&!$([CX3]=[OX1,SX1,NX2])])[NX3]([#1,*])[#1,*]',
	'redux': True,
},
{	
	'name': 'Aminals',
	'description': 'R1, R2, R3, R4, R5, R6 = H, alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([NX3]([#1,*])[#1,*])[NX3]([#1,*])[#1,*]',
	'redux': True,
},
{	
	'name': 'Thiohemiaminals',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([SX2][#1,#6&!$([CX3]=[OX1,SX1,NX2])&!$(C#N)])[NX3]([#1,*])[#1,*]',
	'redux': True,
},
{	
	'name': 'Thioacetals / Thioketals',
	'description': 'R1, R2 = H, alkyl, aryl; R3, R4 = alkyl, aryl',
	'SMARTS': '[CH2,$([CH1][#6]),$([CX4]([#6])[#6])]([SX2][#6&!$([CX3]=[OX1,SX1,NX2])&!$(C#N)])[SX2,OX2][#6&!$([CX3]=[OX1,SX1,NX2])&!$(C#N)]',
	'redux': True,
},
{	
	'name': 'Enamines',
	'description': 'R1, R2, R3, R4, R5 = H, acyl, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[CX3]([#1,#6])[NX3&!$([NX3](=[OX1])=[OX1])&!$([NX3+](=[OX1])[O-])]([#1,*])[#1,*]',
	'redux': True,
},
{	
	'name': 'Enols',
	'description': 'R1, R2, R3 = H, acyl, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[CX3]([OH1])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Enolethers',
	'description': 'R1, R2, R3 = H, acyl, alkyl, aryl; R4 = alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[CX3]([OX2][#6&!$([CX3]([OX2])=[OX1,SX1,NX2,C])])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Hydroxy compounds: alcohols or phenols',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX4]([OH1])[#8,#16,#7,#15])&!$([CX3]([OH1])=[OX1,SX1,NH2,C])][OH1]',
	'redux': True,
},
{	
	'name': 'Alcohols',
	'description': 'R = alkyl',
	'SMARTS': '[$([CX4][OH1])&!$(C([OH1])[#7,#8,#15,#16,F,Cl,Br,I])]',
	'redux': True,
},
{	
	'name': 'Primary alcohols',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[$([#6][CH2][OH1]),$([CH3][OH1])]',
	'redux': True,
},
{	
	'name': 'Secondary alcohols',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[OH1][CX4H1]([#6])[#6]',
	'redux': True,
},
{	
	'name': 'Tertiary alcohols',
	'description': 'R1, R2, R3 = alkyl, aryl',
	'SMARTS': '[CX4]([OH1])([#6])([#6])[#6]',
	'redux': True,
},
{	
	'name': '1,2-Diols',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[CX4]([OH1])([#1,#6])([#1,#6])[CX4]([OH1])([#1,#6])([#1,#6])',
	'redux': True,
},
{	
	'name': '1,2-Aminoalcohols',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[CX4]([#1,#6])([#1,#6])([OH1])[CX4]([#1,#6])([#1,#6])[NH2,$([NH1]([CX4][CX4][OH1])[#6&!$([CX3]=[OX1,SX1,NX2,C])])]',
	'redux': True,
},
{	
	'name': 'Phenols',
	'description': 'Any aromatic or heteroaromatic ring with an OH substituent',
	'SMARTS': 'c[OH1]',
	'redux': True,
},
{	
	'name': 'Diphenols',
	'description': 'Any aromatic or heteroaromatic ring with two OH substituents in ortho-, meta- or para-positions to each other',
	'SMARTS': '[$(c([OH1])c([OH1])),$(c([OH1])[aR1]c([OH1])),$(c1([OH1])aac([OH1])aa1)]',
	'redux': True,
},
{	
	'name': 'Enediols',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([OH1])=[CX3]([#1,#6])([OH1])',
	'redux': True,
},
{	
	'name': 'Ethers',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([OX2])=[OX1,SX1,NX2,C])][OX2][#6&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([OX2])=[OX1,SX1,NX2,C])]',
	'redux': True,
},
{	
	'name': 'Dialkylethers',
	'description': 'R1, R2 = alkyl',
	'SMARTS': '[CX4&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))][OX2][CX4&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))]',
	'redux': True,
},
{	
	'name': 'Alkylarylethers',
	'description': 'R1 = alkyl; R2 = aryl',
	'SMARTS': '[CX4&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))][OX2]c',
	'redux': True,
},
{	
	'name': 'Diarylethers',
	'description': 'R1, R2 = aryl',
	'SMARTS': 'c[OX2&!$([OX2r3])]c',
	'redux': True,
},
{	
	'name': 'Thioethers',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([SX2])=[OX1,SX1,NX2,C])&!$([#6][SX2]C#N)][SX2][#6&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([SX2])=[OX1,SX1,NX2,C])&!$([#6][SX2]C#N)]',
	'redux': True,
},
{	
	'name': 'Disulfides',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[*][#16X2][#16X2][*]',
	'redux': True,
},
{	
	'name': 'Peroxides',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1,NX2])][OX2][OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Hydroperoxides',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#1,#6][OX2][OX2H1]',
	'redux': True,
},
{	
	'name': 'Hydrazine derivatives',
	'description': 'R1, R2, R3, R4 = H, acyl, alkyl, aryl',
	'SMARTS': '[$([NX4]([#1,#6])([#1,#6])([#1,#6])),$([NX3]([#1,#6])([#1,#6])),$([NX2]=[#1,#6])][$([NX4]([#1,#6])([#1,#6])([#1,#6])),$([NX3]([#1,#6])([#1,#6])),$([NX2]=[#1,#6])]',
	'redux': True,
},
{	
	'name': 'Hydroxylamines',
	'description': 'R1, R2, R3 = H, alkyl, aryl',
	'SMARTS': '[N&!$([NH1]([CX3]=[OX1,SX1,NX2])[OH1])&!$(N([OH1])([CX3]=[OX1,SX1,NX2])([CX3]=[OX1,SX1,NX2]))]([#1,#6])([#1,#6])[OH1]',
	'redux': True,
},
{	
	'name': 'Amines',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl',
	'SMARTS': '[NX3]([#1,#6&!$([CX3]=[O,S,N,C])&!$([#6]([NX3])[N,O,S,P])&!$(C([NX3])#N)])([#1,#6&!$([CX3]=[O,S,N,C])&!$([#6]([NX3])[N,O,S,P])&!$(C([NX3])#N)])[#1,#6&!$([CX3]=[O,S,N,C])&!$([#6]([NX3])[N,O,S,P])&!$(C([NX3])#N)]',
	'redux': True,
},
{	
	'name': 'Primary aliphatic amines',
	'description': 'R = alkyl',
	'SMARTS': '[NX3H2][CX4&!$([CX4]([NH2])[O,N,S,P])]',
	'redux': True,
},
{	
	'name': 'Primary aromatic amines',
	'description': 'R = aryl',
	'SMARTS': '[NX3H2]c',
	'redux': True,
},
{	
	'name': 'Secondary aliphatic amines',
	'description': 'R1, R2 = alkyl',
	'SMARTS': '[NX3H1]([CX4&!$([CX4]([NX3H1])[O,S,N,P])])[CX4&!$([CX4]([NX3H1])[O,S,N,P])]',
	'redux': True,
},
{	
	'name': 'Secondary mixed amines (aryl alkyl)',
	'description': 'R1 = alkyl; R2 = aryl',
	'SMARTS': 'c[NX3H1][CX4&!$([CX4]([NX3H1])[O,N,S,P])]',
	'redux': True,
},
{	
	'name': 'Secondary aromatic amines',
	'description': 'R1, R2 = aryl',
	'SMARTS': '[NX3H1](c)c',
	'redux': True,
},
{	
	'name': 'Tertiary aliphatic amines',
	'description': 'R1, R2, R3 = alkyl',
	'SMARTS': '[NX3]([CX4&!$([CX4]([NX3])[O,S,N,P])])([CX4&!$([CX4]([NX3])[O,S,N,P])])[CX4&!$([CX4]([NX3])[O,S,N,P])]',
	'redux': True,
},
{	
	'name': 'Tertiary mixed amines (aryl alkyl)',
	'description': 'R1 = alkyl; R2 = aryl; R3 = alkyl, aryl',
	'SMARTS': 'c[NX3]([CX4&!$([CX4]([NX3])[O,N,S,P])])[c,CX4&!$([CX4]([NX3]c)[O,S,N,P])]',
	'redux': True,
},
{	
	'name': 'Tertiary aromatic amines',
	'description': 'R1, R2, R3 = aryl',
	'SMARTS': '[NX3](c)(c)c',
	'redux': True,
},
{	
	'name': 'Quaternary ammonium salts',
	'description': 'R1, R2, R3, R4 = alkyl, aryl',
	'SMARTS': '[NX4,NX4+]([#6])([#6])([#6])[#6]',
	'redux': True,
},
{	
	'name': 'N-Oxides',
	'description': 'R1, R2, R3 = alkyl, aryl',
	'SMARTS': '[$([OX1-,OH1][#7X4+]([*])([*])([*])),$([OX1]=[#7X4v5]([*])([*])([*])),$([OX1]=[Nv5X3&!$([Nv5X3](=[OX1])=[OX1])](=[*])([*])),$([OX1-,OH1][#7X3+R](~[R])(~[R])),$([OX1]=[#7v5X3R](~[R])(~[R]))]',
	'redux': True,
},
{	
	'name': 'Halogen derivatives (alkyl, alkenyl, aryl)',
	'description': 'R = alkyl, alkenyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[CX4,CX3,c;!$([#6]=[OX1])][F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Alkyl fluorides',
	'description': 'R = alkyl',
	'SMARTS': '[CX4]F AND NOT [$([CX4](F)[CX3]=[CX3]),$([CX4](F)[CX2]#[CX2]),$([CX4](F)a)]',
	'redux': True,
},
{	
	'name': 'Alkyl chlorides',
	'description': 'R = alkyl',
	'SMARTS': '[CX4]Cl AND NOT [$([CX4](Cl)[CX3]=[CX3]),$([CX4](Cl)[CX2]#[CX2]),$([CX4](Cl)a)]',
	'redux': True,
},
{	
	'name': 'Alkyl bromides',
	'description': 'R = alkyl',
	'SMARTS': '[CX4]Br AND NOT [$([CX4](Br)[CX3]=[CX3]),$([CX4](Br)[CX2]#[CX2]),$([CX4](Br)a)]',
	'redux': True,
},
{	
	'name': 'Alkyl iodides',
	'description': 'R = alkyl',
	'SMARTS': '[CX4]I AND NOT [$([CX4](I)[CX3]=[CX3]),$([CX4](I)[CX2]#[CX2]),$([CX4](I)a)]',
	'redux': True,
},
{	
	'name': 'Aryl fluorides',
	'description': 'R = aryl',
	'SMARTS': 'aF',
	'redux': True,
},
{	
	'name': 'Aryl chlorides',
	'description': 'R = aryl',
	'SMARTS': 'aCl',
	'redux': True,
},
{	
	'name': 'Aryl bromides',
	'description': 'R = aryl',
	'SMARTS': 'aBr',
	'redux': True,
},
{	
	'name': 'Aryl iodides',
	'description': 'R = aryl',
	'SMARTS': 'aI',
	'redux': True,
},
{	
	'name': 'Organometallic compounds',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#32!#33!#34!#35!#36!#51!#52!#53!#54!#84!#85!#86]~[#6]',
	'redux': True,
},
{	
	'name': 'Organolithium compounds',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[Li][#6]',
	'redux': True,
},
{	
	'name': 'Organomagnesium compounds',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[Mg][#6]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid derivatives',
	'description': 'R = H, alkyl, aryl; Y = any hetero atom',
	'SMARTS': '[#1,#6][CX3](=[OX1])[!#1!#6!$([SH1,SX2])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acids',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[OX2H1]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid salts',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[OX1-]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid esters',
	'description': 'R1 = H, alkyl, aryl; R2 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Lactones',
	'description': 'A cyclic ester (any ring size)',
	'SMARTS': '[#6R][CX3R](=[OX1])[OX2R][#6R&!$([CX3R]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid amides',
	'description': 'R1, R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX1])])[#1,#6&!$([CX3]=[OX1,SX1,NX1])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid primary amides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[NX3H2]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid secondary amides',
	'description': 'R1 = H, alkyl, aryl; R2 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[NX3H1][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid tertiary amides',
	'description': 'R1 = H, alkyl, aryl; R2, R3 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[NX3]([#6&!$([CX3]=[OX1,SX1,NX2])])[#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Lactams',
	'description': 'R = H, alkyl, aryl; a cyclic amide (any ring size)',
	'SMARTS': '[#6R][CX3R](=[OX1])[NX3R]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#6R&!$([CX3R]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid hydrazides',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[$([#6,#1][#6X3](=[OX1])[#7X3&!$([#7X3]([*]=[OX1,SX1])[#7X3]([*]=[OX1,SX1]))][#7X3]),$(C1(=O)[NX3][NX3]C(=O)cc1)]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid azides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[$([NX2]=[NX2+]=[NX1-]),$([NX2]=[NX2+]=N),$([NX2-]-[NX2+]#[NX1])]',
	'redux': True,
},
{	
	'name': 'Hydroxamic acids',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[NX3][OH1]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid amidines',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][NX2]=[CX3]([#1,#6])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,N])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid amidrazones',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][NX2]=[CX3]([#1,#6])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#1,#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Nitriles',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6&!$([CX3]=[OX1,SX1])][CX2]#[NX1]',
	'redux': True,
},
{	
	'name': 'Acyl halides',
	'description': 'R = H, alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#1,#6][CX3](=[OX1])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Acyl fluorides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])F',
	'redux': True,
},
{	
	'name': 'Acyl chlorides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])Cl',
	'redux': True,
},
{	
	'name': 'Acyl bromides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])Br',
	'redux': True,
},
{	
	'name': 'Acyl iodides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])I',
	'redux': True,
},
{	
	'name': 'Acyl cyanides',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[OX1])[CX2]#[NX1]',
	'redux': True,
},
{	
	'name': 'Imido esters',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[NX2][#1,#6])[OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])]',
	'redux': True,
},
{	
	'name': 'Imidoyl halides',
	'description': 'R1, R2 = H, alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#1,#6][CX3](=[NX2][#1,#6])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Thiocarboxylic acid derivatives',
	'description': 'R = H, alkyl, aryl; Y = any hetero atom',
	'SMARTS': '[$([#1,#6][CX3](=[SX1])[!#1!#6]),$([#1,#6][CX3](=[OX1])[SX2])]',
	'redux': True,
},
{	
	'name': 'Thiocarboxylic acids',
	'description': 'R = H, alkyl, aryl; Y = OH, SH',
	'SMARTS': '[$([#1,#6][CX3](=[SX1])[OH1,SH1]),$([#1,#6][CX3](=[OX1])[SH1])]',
	'redux': True,
},
{	
	'name': 'Thiocarboxylic acid esters',
	'description': 'R1 = H, alkyl, aryl; R2 = aklyl, aryl',
	'SMARTS': '[$([#1,#6][CX3](=[SX1])[OX2,SX2][#6&!$([CX3]=[OX1,SX1,NX2])]),$([#1,#6][CX3](=[OX1])[SX2][#6&!$([CX3]=[OX1,SX1,NX2])])]',
	'redux': True,
},
{	
	'name': 'Thiolactones',
	'description': 'A cyclic thioester (any ring size)',
	'SMARTS': '[$([#6R][CX3R](=[SX1])[OX2R,SX2R][#6R&!$([CX3R]=[OX1,SX1,NX2])]),$([#6R][CX3R](=[OX1])[SX2R][#6R&!$([CX3R]=[OX1,SX1,NX2])])]',
	'redux': True,
},
{	
	'name': 'Thiocarboxylic acid amides',
	'description': 'R1, R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[SX1])[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX1])])[#1,#6&!$([CX3]=[OX1,SX1,NX1])]',
	'redux': True,
},
{	
	'name': 'Thiolactams',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[#6R][CX3R](=[SX1])[NX3R]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])[#6R&!$([CX3R]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Imidothioesters',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX3](=[NX2][#1,#6])[SX2][C&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Oxohetarenes',
	'description': 'R = H, alkyl, aryl; any heteroaromatic compound with a C=O structure',
	'SMARTS': '[$(c(=[OX1])[$(n[#1,#6]),o,s]),$(c(=[OX1])a[$(n[#1,#6]),o,s]),$(c(=[OX1])aa[$(n[#1,#6]),o,s]),$(c(=[OX1])aaa[$(n[#1,#6]),o,s])]',
	'redux': True,
},
{	
	'name': 'Thioxohetarenes',
	'description': 'R = H, alkyl, aryl; any heteroaromatic compound with a C=S structure',
	'SMARTS': '[$(c(=[SX1])[$(n[#1,#6]),o,s]),$(c(=[SX1])a[$(n[#1,#6]),o,s]),$(c(=[SX1])aa[$(n[#1,#6]),o,s]),$(c(=[SX1])aaa[$(n[#1,#6]),o,s])]',
	'redux': True,
},
{	
	'name': 'Iminohetarenes',
	'description': 'R1, R2 = H, alkyl, aryl; ; any heteroaromatic compound with a C=N structure',
	'SMARTS': '[$(c(=N)[$(n[#1,#6]),o,s]),$(c(=N)a[$(n[#1,#6]),o,s]),$(c(=N)aa[$(n[#1,#6]),o,s]),$(c(=N)aaa[$(n[#1,#6]),o,s])]',
	'redux': True,
},
{	
	'name': 'Orthocarboxylic acid derivatives',
	'description': 'R = H, alkyl, aryl; Y = OH, alkoxy, aryloxy, (substituted) amino, etc.',
	'SMARTS': '[CX4&!$([CX4]([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I]))]([#1,#6&!$([CX3]=[OX1,SX1,NX2])])([!#1!#6])([!#1!#6])[!#1!#6]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid orthoesters',
	'description': 'R1 = H, alkyl, aryl; R2, R3, R4 = alkyl, aryl',
	'SMARTS': '[#1,#6][CX4]([OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])])([OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])])[OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid amide acetals',
	'description': 'R1, R2, R3 = H, alkyl, aryl; R4, R5 = alkyl, aryl',
	'SMARTS': '[CH1,$([CX4][#6])]([OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])])([OX2][#6&!$([CX3]=[OX1,SX1,NX2,C])])[NX3]([#1,*])[#1,*]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid anhydrides',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[#6][#6X3](=[OX1])[#8X2][#6X3](=[OX1])[#6]',
	'redux': True,
},
{	
	'name': 'Carboxylic acid imides',
	'description': 'R1, R2, R3 = H, alkyl, aryl',
	'SMARTS': '[CX3](=[OX1])[NX3][CX3](=[OX1])',
	'redux': True,
},
{	
	'name': 'Carboxylic acid unsubstituted imides',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX3](=[OX1])[NX3H1][CX3](=[OX1])',
	'redux': True,
},
{	
	'name': 'Carboxylic acid substituted imides',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = any atom except H',
	'SMARTS': '[CX3](=[OX1])[NX3H0][CX3](=[OX1])',
	'redux': True,
},

{
	'name': 'CO2 derivative (general)',
	'description': 'Any carbon with 4 valences to hetero atoms; A = any atom except H and carbon',
	'SMARTS': '[$([Cv4X4]([!#1!#6])([!#1!#6])([!#1!#6])[!#1!#6]),$([CX3](=[!#1!#6])([!#1!#6])[!#1!#6]),$([CX2](=[!#1!#6])=[!#1!#6]),$([CX2](#[!#1!#6])[!#1!#6])]',
	'redux': True,
},
{	
	'name': 'Carbonic acid derivatives',
	'description': 'R1, R2 = any hetero atom',
	'SMARTS': '[$([CX3](=[OX1])([!#1!#6])[!#1!#6])&!$([CX3](=[OX1])([#7,#16X2,#15,#5])[#8X2,#16X2,#7,#15,#5])]',
	'redux': True,
},
{	
	'name': 'Carbonic acid monoesters',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[OH1][CX3](=[OX1])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Carbonic acid diesters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1,NX2])][OX2][CX3](=[OX1])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Carbonic acid ester halides',
	'description': 'R = alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#6][OX2][CX3](=[OX1])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Thiocarbonic acid derivatives',
	'description': 'R1, R2 = any hetero atom',
	'SMARTS': '[$([CX3](=[SX1])([!#1!#6])[!#1!#6])&!$([CX3](=[SX1])([#7,#16X2,#15,#5])[#8X2,#16X2,#7,#15,#5]),$([CX3](=[OX1])([!#1!#6])[SX2])&!$([CX3](=[OX1])([#7,#15,#5])[#8X2,#16X2,#7,#15,#5])]',
	'redux': True,
},
{	
	'name': 'Thiocarbonic acid monoesters',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[$([OH1,SH1][CX3](=[SX1])[OX2,SX2][#6&!$([CX3]=[OX1,SX1,NX2])]),$([OH1][CX3](=[OX1])[SX2][#6&!$([CX3]=[OX1,SX1,NX2])]),$([SH1][CX3](=[SX1,OX1])[OX2,SX2][#6&!$([CX3]=[OX1,SX1,NX2])])]',
	'redux': True,
},
{	
	'name': 'Thiocarbonic acid diesters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[$([#6&!$([CX3]=[OX1,SX1,NX2])][OX2,SX2][CX3](=[SX1])[OX2,SX2][#6&!$([CX3]=[OX1,SX1,NX2])]),$([#6&!$([CX3]=[OX1,SX1,NX2])][SX2][CX3](=[OX1])[OX2,SX2][#6&!$([CX3]=[OX1,SX1,NX2])])]',
	'redux': True,
},
{	
	'name': 'Thiocarbonic acid ester halides',
	'description': 'R = alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[$([#6&!$([CX3]=[OX1,SX1,NX2])][OX2,SX2][CX3](=[SX1])[F,Cl,Br,I]),$([#6&!$([CX3]=[OX1,SX1,NX2])][SX2][CX3](=[OX1])[F,Cl,Br,I])]',
	'redux': True,
},
{	
	'name': 'Carbamic acid derivatives',
	'description': 'R1, R2 = H, alkyl, aryl; Y = OH, alkoxy, aryloxy, halogen',
	'SMARTS': '[NX3]([#1,*])([#1,*])[CX3](=[OX1])[!#1!#6!#16!$([Pv5])]',
	'redux': True,
},
{	
	'name': 'Carbamic acid',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[OH1][CX3](=[OX1])[NX3]',
	'redux': True,
},
{	
	'name': 'Carbamic acid esters (urethanes)',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[#6][OX2][CX3](=[OX1])[NH2,NH1,NX3]',
	'redux': True,
},
{	
	'name': 'Carbamic acid halides',
	'description': 'R1, R2 = H, alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[F,Cl,Br,I][CX3](=[OX1])[NX3]',
	'redux': True,
},
{	
	'name': 'Thiocarbamic acid derivatives',
	'description': 'R1, R2 = H, alkyl, aryl; Y = OH, alkoxy, aryloxy, halogen',
	'SMARTS': '[$([NX3&!$([NX3]([CX3]=[SX1])[CX3]=[OX1,SX1,NX2])][CX3](=[SX1])[!#1!#6!#7!$([Pv5])]),$([NX3&!$([NX3]([CX3]=[OX1])[CX3]=[OX1,SX1,NX2])][CX3](=[OX1])[SH1,SX2])]',
	'redux': True,
},
{	
	'name': 'Thiocarbamic acids',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[$([OH1,SH1][CX3](=[SX1])[NX3]),$([SH1][CX3](=[OX1])[NX3])]',
	'redux': True,
},
{	
	'name': 'Thiocarbamic acid esters',
	'description': 'R1, R2 = H, alkyl, aryl; R3 = alkyl, aryl',
	'SMARTS': '[$([#6][OX2,SX2][CX3](=[SX1])[NX3&!$([NX3]([CX3]=[OX1,SX1])[CX3]=[OX1,SX1])]),$([#6][SX2][CX3](=[OX1])[NX3&!$([NX3]([CX3]=[OX1,SX1])[CX3]=[OX1,SX1])])]',
	'redux': True,
},
{	
	'name': 'Thiocarbamic acid halides',
	'description': 'R1, R2 = H, alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[F,Cl,Br,I][CX3](=[SX1])[NX3]',
	'redux': True,
},
{	
	'name': 'Ureas',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[#7&!$([#7]=[#7])&!$([#7][#7H2])&!$([#7][#7H1][#6])&!$([#7][#7]([#1,#6])[#1,#6])&!$([#7][#7]=[#6])][CX3](=[OX1])[#7&!$([#7]=[#7])&!$([#7][#7H2])&!$([#7][#7H1][#6])&!$([#7][#7]([#1,#6])[#1,#6])&!$([#7][#7]=[#6])]',
	'redux': True,
},
{	
	'name': 'Isoureas',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[NX3&!$([NX3]([CX3]=[OX1,SX1,NX2])([CX3]=[OX1,SX1,NX2]))]([#1,#6])([#1,#6])[CX3]([OX2])=[NX2][#1,#6]',
	'redux': True,
},
{	
	'name': 'Thioureas',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[#7&!$([#7]=[#7])&!$([#7][#7H2])&!$([#7][#7H1][#6])&!$([#7][#7]([#1,#6])[#1,#6])&!$([#7][#7]=[#6])][CX3](=S)[#7&!$([#7]=[#7])&!$([#7][#7H2])&!$([#7][#7H1][#6])&!$([#7][#7]([#1,#6])[#1,#6])&!$([#7][#7]=[#6])]',
	'redux': True,
},
{	
	'name': 'Isothioureas',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[NX3&!$([NX3]([CX3]=[OX1,SX1,NX2])([CX3]=[OX1,SX1,NX2]))]([#1,#6])([#1,#6])[CX3]([S&!$([SX2]([CX3]=N)C=C)])=[NX2][#1,#6]',
	'redux': True,
},
{	
	'name': 'Guanidines',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': 'N[CX3](=[NX2][#1,#6])N',
	'redux': True,
},
{	
	'name': 'Semicarbazides',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[$([NX3][CX3](=[OX1])[NX3&!$([NX3]([CX3]=[OX1,SX1])[NX3][CX3]=[OX1,SX1])][NX3])&!$([NX3][CX3](=[OX1])[NX3]N=C)]',
	'redux': True,
},
{	
	'name': 'Thiosemicarbazides',
	'description': 'R1, R2, R3, R4, R5 = H, alkyl, aryl',
	'SMARTS': '[$([NX3][CX3](=[SX1])[NX3&!$([NX3]([CX3]=[OX1,SX1])[NX3][CX3]=[OX1,SX1])][NX3])&!$([NX3][CX3](=[SX1])[NX3]N=C)]',
	'redux': True,
},
{	
	'name': 'Azides',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][$([NX2]=[NX2+]=[NX1-]),$(N=[NX2+]=N),$([NX2-]-[NX2+]#[NX1]),$([NX2]=[NX2]#[NX1])]',
	'redux': True,
},
{	
	'name': 'Azo compounds',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][NX2]=[NX2][#6]',
	'redux': True,
},
{	
	'name': 'Diazonium salts',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][NX2+]#[NX1]',
	'redux': True,
},
{	
	'name': 'Isonitriles',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][$([NX2+]#[CX1-]),$([NX2]#[CX1])]',
	'redux': True,
},
{	
	'name': 'Cyanates',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][OX2][CX2]#[NX1]',
	'redux': True,
},
{	
	'name': 'Isocyanates',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][NX2]=[CX2]=[OX1]',
	'redux': True,
},
{	
	'name': 'Thiocyanates',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][SX2][CX2]#[NX1]',
	'redux': True,
},
{	
	'name': 'Isothiocyanates',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][NX2]=[CX2]=[SX1]',
	'redux': True,
},
{	
	'name': 'Carbodiimides',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[#1,#6&!$([CX3]=[OX1,SX1])][NX2]=[CX2]=[NX2][#1,#6&!$([CX3]=[OX1,SX1])]',
	'redux': True,
},
{	
	'name': 'Nitroso compounds',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6,#7][NX2]=[OX1]',
	'redux': True,
},
{	
	'name': 'Nitro compounds',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-])]',
	'redux': True,
},
{	
	'name': 'Nitrites',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][OX2][NX2]=[OX1]',
	'redux': True,
},
{	
	'name': 'Nitrates',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1])][OX2][$([NX3](=[OX1])=[OX1]),$([NX3+](=[OX1])[O-])]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid derivatives',
	'description': 'Y, Z = any hetero atom',
	'SMARTS': '[!#1!#6][Sv6X4](=[OX1])(=[OX1])[!#1!#6]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid',
	'description': '',
	'SMARTS': '[OH1,O-][Sv6X4](=[OX1])(=[OX1])[OX2,O-]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid monoesters',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][OX2][Sv6X4](=[OX1])(=[OX1])[OH1,O-]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid diesters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][OX2][Sv6X4](=[OX1])(=[OX1])[OX2][#6]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid amide esters',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#6][OX2][Sv6X4](=[OX1])(=[OX1])N',
	'redux': True,
},
{	
	'name': 'Sulfuric acid amides',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[OH1,O-][Sv6X4](=[OX1])(=[OX1])N',
	'redux': True,
},
{	
	'name': 'Sulfuric acid diamides',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': 'N[Sv6X4](=[OX1])(=[OX1])N',
	'redux': True,
},
{	
	'name': 'Sulfuryl halides',
	'description': 'X = F, Cl, Br, I; Y = any hetero atom',
	'SMARTS': '[!#1!#6][Sv6X4](=[OX1])(=[OX1])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Sulfonic acid derivatives',
	'description': 'R = alkyl, aryl; Y = any hetero atom',
	'SMARTS': '[#6][Sv6X4](=[OX1])(=[OX1])[!#1!#6]',
	'redux': True,
},
{	
	'name': 'Sulfonic acids',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][Sv6X4](=[OX1])(=[OX1])[OH1,O-]',
	'redux': True,
},
{	
	'name': 'Sulfonic acid esters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][Sv6X4](=[OX1])(=[OX1])[OX2][#6]',
	'redux': True,
},
{	
	'name': 'Sulfonamides',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#6][Sv6X4](=[OX1])(=[OX1])N',
	'redux': True,
},
{	
	'name': 'Sulfonyl halides',
	'description': 'R = alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#6][Sv6X4](=[OX1])(=[OX1])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Sulfones',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][Sv6X4](=[OX1,SX1,NX2])(=[OX1,SX1,NX2])[#6]',
	'redux': True,
},
{	
	'name': 'Sulfoxides',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][Sv4X3](=[OX1])[#6]',
	'redux': True,
},
{	
	'name': 'Sulfinic acid derivatives',
	'description': 'R = alkyl, aryl; Y = any hetero atom',
	'SMARTS': '[#6][Sv4X3](=[OX1])[!#1!#6]',
	'redux': True,
},
{	
	'name': 'Sulfinic acids',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6][Sv4X3](=[OX1])[OH1,O-]',
	'redux': True,
},
{	
	'name': 'Sulfinic acid esters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6][Sv4X3](=[OX1])[OX2][#6]',
	'redux': True,
},
{	
	'name': 'Sulfinic acid halides',
	'description': 'R = alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#6][Sv4X3](=[OX1])[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Sulfinic acid amides',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#6][Sv4X3](=[OX1])[NX3&!$([NX3][CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Sulfenic acid derivatives',
	'description': 'R = alkyl, aryl; Y = any hetero atom',
	'SMARTS': '[#6][Sv2X2][!#6!#1!S]',
	'redux': True,
},
{	
	'name': 'Sulfenic acids',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[SX1,OX1,NX2])][Sv2X2][OH1]',
	'redux': True,
},
{	
	'name': 'Sulfenic acid esters',
	'description': 'R1, R2 = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[SX1,OX1,NX2])][Sv2X2][OX2][#6]',
	'redux': True,
},
{	
	'name': 'Sulfenic acid halides',
	'description': 'R = alkyl, aryl; X = F, Cl, Br, I',
	'SMARTS': '[#6&!$([CX3]=[SX1,OX1,NX2])][Sv2X2][F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Sulfenic acid amides',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[SX1,OX1,NX2])][Sv2X2][NX3]([#6&!$([CX3]=[OX1,SX1,NX2])])[#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Thiols',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[#6&!$([CX3]=[OX1,SX1,N,C])&!$(C#N)][Sv2H1]',
	'redux': True,
},
{	
	'name': 'Alkylthiols',
	'description': 'R = alkyl',
	'SMARTS': '[C&!$([CX3]=[OX1,SX1,N,C])&!$(C#N)][Sv2H1]',
	'redux': True,
},
{	
	'name': 'Thiophenols',
	'description': 'R = aryl',
	'SMARTS': 'c[Sv2H1]',
	'redux': True,
},
{	
	'name': 'Phosphoric acid derivatives',
	'description': 'X, Y, Z = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([H,!#6!$([SX2])])([H,!#6!$([SX2])])[H,!#6!$([SX2])]',
	'redux': True,
},
{	
	'name': 'Phosphoric acids',
	'description': '',
	'SMARTS': '[OX1]=[Pv5X4]([OH1])([O])[O]',
	'redux': True,
},
{	
	'name': 'Phosphoric acid esters',
	'description': 'R = alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([OX2][#6])([!#6!$([SX2])])[!#6!$([SX2])]',
	'redux': True,
},
{	
	'name': 'Phosphoric acid halides',
	'description': 'X = F, Cl, Br, I; Y, Z = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([F,Cl,Br,I])([!#6!$([SX2])])[!#6!$([SX2])]',
	'redux': True,
},
{	
	'name': 'Phosphoric acid amides',
	'description': 'R1, R2 = H, alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([NX3])([!#6!$([SX2])])[!#6!$([SX2])]',
	'redux': True,
},
{	
	'name': 'Thiophosphoric acid derivatives',
	'description': 'X, Y, Z = any O, N, Hal residue',
	'SMARTS': '[$([SX1]=[Pv5X4]([H,!#6])([H,!#6])[H,!#6]),$([OX1]=[Pv5X4]([SX2])([H,!#6])[H,!#6])]',
	'redux': True,
},
{	
	'name': 'Thiophosphoric acids',
	'description': '',
	'SMARTS': '[$([SX1]=[Pv5X4]([OH1,SH1])([O,S])[O,S]),$([OX1]=[Pv5X4]([SH1])([O,S])[O,S])]',
	'redux': True,
},
{	
	'name': 'Thiophosphoric acid esters',
	'description': 'R = alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[$([SX1]=[Pv5X4]([OX2,SX2][#6])([!#6])[!#6]),$([OX1]=[Pv5X4]([SX2][#6])([!#6])[!#6])]',
	'redux': True,
},
{	
	'name': 'Thiophosphoric acid halides',
	'description': 'X = F, Cl, Br, I; Y, Z = any O, N, Hal residue',
	'SMARTS': '[$([SX1]=[Pv5X4]([F,Cl,Br,I])([!#6])[!#6]),$([OX1]=[Pv5X4]([F,Cl,Br,I])([SX2,SH1])[!#6])]',
	'redux': True,
},
{	
	'name': 'Thiophosphoric acid amides',
	'description': 'R1, R2 = H, alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[$([SX1]=[Pv5X4]([NX3])([!#6])[!#6]),$([OX1]=[Pv5X4]([NX3])([SX2,SH1])[!#6])]',
	'redux': True,
},
{	
	'name': 'Phosphonic acid derivatives',
	'description': 'R = alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([#6])([!#6])[!#6]',
	'redux': True,
},
{	
	'name': 'Phosphonic acids',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[OX1]=[Pv5X4]([#6])([OH1])O',
	'redux': True,
},
{	
	'name': 'Phosphonic acid esters',
	'description': 'R1, R2 = alkyl, aryl; Y = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([#6])([OX2][#6])[!#6]',
	'redux': True,
},
{	
	'name': 'Phosphines',
	'description': 'R1, R2, R3 = alkyl, aryl',
	'SMARTS': '[Pv3X3]([#6])([#6])[#6]',
	'redux': True,
},
{	
	'name': 'Phosphinoxides',
	'description': 'R1, R2, R3 = alkyl, aryl',
	'SMARTS': '[Pv5X4](=[OX1])([#6])([#6])[#6]',
	'redux': True,
},
{	
	'name': 'Boronic acid derivatives',
	'description': 'R = alkyl, aryl; Y, Z = any O, N, Hal residue',
	'SMARTS': '[Bv3X3]([#6])([!#6])[!#6]',
	'redux': True,
},
{	
	'name': 'Boronic acids',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[BX3]([#6])([OH1,SH1])[O,S]',
	'redux': True,
},
{	
	'name': 'Boronic acid esters',
	'description': 'R1, R2 = alkyl, aryl; Y = any O, N, Hal residue',
	'SMARTS': '[BX3]([#6])([OX2][#6])[H,!#6]',
	'redux': True,
},
{	
	'name': 'Alkynes',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Aromatic compounds',
	'description': 'Any aromatic carbocyclic or heterocyclic structure',
	'SMARTS': '[a!R0]',
	'redux': True,
},
{	
	'name': 'Arenes',
	'description': 'Any compounds containing at least one benzen ring',
	'SMARTS': 'c1ccccc1',
	'redux': True,
},
{	
	'name': 'Heterocyclic compounds',
	'description': 'Any cyclic structure with at least one non-carbon atom incorporated',
	'SMARTS': '[*R]~[*R&!#6]~[*R]',
	'redux': True,
},
{	
	'name': 'Aromatic heterocyclic compounds',
	'description': 'Any aromatic cyclic structure with at least one non-carbon atom incorporated',
	'SMARTS': '[a!#6!R0]',
	'redux': True,
},
{	
	'name': 'alpha-Aminoacids',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[OH1][CX3](=[OX1])[CX4]([#1,#6])([#1,#6])[NH2,$([NH1]([#6])[#6])&!$([NH1](CC(=O)[OH1])[CX3]=[OX1,SX1])]',
	'redux': True,
},
{	
	'name': 'alpha-Hydroxyacids',
	'description': 'R = H, alkyl, aryl',
	'SMARTS': '[OH1][CX3](=[OX1])[CX4]([OH1])[#1,#6]',
	'redux': True,
},
{	
	'name': '1,2-Diamines',
	'description': 'R1, R2, R3, R4 = H , any carbon atom; R5, R6, R7, R8 = H, any carbon atom except Csp2 in amide and thioamide groups',
	'SMARTS': '[NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])[CX4][CX4][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])',
	'redux': True,
},
{	
	'name': '1,2-Dithiols',
	'description': 'R1, R2, R3, R4 = H , any carbon atom',
	'SMARTS': '[CX4]([SH1])([#1,#6])([#1,#6])[CX4]([SH1])([#1,#6])([#1,#6])',
	'redux': True,
},
{	
	'name': '1,2-Aminothiols',
	'description': 'R1, R2, R3, R4 = H , any carbon atom; R5, R6 = H, any carbon atom except Csp2 in amide and thioamide groups',
	'SMARTS': '[SH1][CX4][CX4][NX3]([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])([#1,#6&!$([CX3]=[OX1,SX1,NX2,C])])',
	'redux': True,
},
{	
	'name': 'alpha-Oxoacids',
	'description': 'R = H, any carbon atom',
	'SMARTS': '[OH1][CX3](=[OX1])[CX3](=[OX1])[#1,#6]',
	'redux': True,
},
{	
	'name': 'alpha, beta-Unsaturated carboxylic acids',
	'description': 'R1, R2, R3 = H, any carbon atom',
	'SMARTS': '[OH1][CX3](=[OX1])[CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Dithiophenols',
	'description': 'Any aromatic or heteroaromatic ring with two SH substituents in ortho-, meta- or para-positions to each other',
	'SMARTS': '[$(c([SH1])c([SH1])),$(c([SH1])[aR1]c([SH1])),$(c1([SH1])aac([SH1])aa1)]',
	'redux': True,
},
{	
	'name': 'Aminophenols',
	'description': 'Any aromatic or heteroaromatic ring with amino group in ortho-, meta- or para-positions to OH group; R1, R2 = H, any carbon atom',
	'SMARTS': '[$(c([OH1])c[NX3]([#1,#6])[#1,#6]),$(c([OH1])[aR1]c[NX3]([#1,#6])[#1,#6]),$(c1([OH1])aac([NX3]([#1,#6])[#1,#6])aa1)]',
	'redux': True,
},
{	
	'name': 'Aminothiophenols',
	'description': 'Any aromatic or heteroaromatic ring with amino group in ortho-, meta- or para-positions to SH group; R1, R2 = H, any carbon atom',
	'SMARTS': '[$(c([SH1])c[NX3]([#1,#6])[#1,#6]),$(c([SH1])[aR1]c[NX3]([#1,#6])[#1,#6]),$(c1([SH1])aac([NX3]([#1,#6])[#1,#6])aa1)]',
	'redux': True,
},
{	
	'name': 'Aryl halides',
	'description': 'R = aryl',
	'SMARTS': 'a[F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Alkyl halides',
	'description': 'R = alkyl',
	'SMARTS': '[CX4][F,Cl,Br,I] AND NOT [$([CX4]([F,Cl,Br,I])[CX3]=[CX3]),$([CX4]([F,Cl,Br,I])[CX2]#[CX2]),$([CX4]([F,Cl,Br,I])a)]',
	'redux': True,
},
{	
	'name': 'Primary amines',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[NX3H2][c,CX4&!$([CX4]([NH2])[O,N,S,P])]',
	'redux': True,
},
{	
	'name': 'Secondary amines',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[NX3H1]([c,CX4&!$([CX4]([NX3H1])[O,S,N,P])])[c,CX4&!$([CX4]([NX3H1])[O,S,N,P])]',
	'redux': True,
},
{	
	'name': 'Tertiary amines',
	'description': 'R = alkyl, aryl',
	'SMARTS': '[NX3]([c,CX4&!$([CX4]([NX3])[O,S,N,P])])([c,CX4&!$([CX4]([NX3])[O,S,N,P])])[c,CX4&!$([CX4]([NX3])[O,S,N,P])]',
	'redux': True,
},
{	
	'name': 'Thiocarbonic acid esters',
	'description': 'R1 = alkyl, aryl; R2 = H, alkyl, aryl',
	'SMARTS': '[$([#6&!$([CX3]=[OX1,SX1,NX2])][OX2,SX2][CX3](=[SX1])[OX2,SX2][#1,#6&!$([CX3]=[OX1,SX1,NX2])]),$([#6&!$([CX3]=[OX1,SX1,NX2])][SX2][CX3](=[OX1])[OX2,SX2][#1,#6&!$([CX3]=[OX1,SX1,NX2])]),$([#6&!$([CX3]=[OX1,SX1,NX2])][OX2,SX2][CX3](=[OX1])[SH1])]',
	'redux': True,
},
{	
	'name': 'Acid anhydrides',
	'description': 'R1, R2 = any atom/ group',
	'SMARTS': '[#6X3](=[OX1])[#8X2][#6X3](=[OX1])',
	'redux': True,
},
{	
	'name': 'Conjugated alkadienes (1,3-alkadienes)',
	'description': 'R1, R2, R3, R4, R5, R6 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Unconjugated alkadienes (1,4-alkadienes)',
	'description': 'R1, R2, R3, R4, R5, R6, R7, R8 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4][CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Vinyl halides',
	'description': 'R1, R2, R3 = any atom/ group; X = F, Cl, Br, I',
	'SMARTS': '[F,Cl,Br,I][CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Gem-trihalides',
	'description': 'R = any atom/ group except halogen atoms; X = F, Cl, Br, I',
	'SMARTS': '[CX4]([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I]) AND NOT [CX4]([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I])([F,Cl,Br,I])',
	'redux': True,
},
{	
	'name': 'Alkynyl halides',
	'description': 'R = any atom/ group; X = F, Cl, Br, I',
	'SMARTS': '[F,Cl,Br,I][CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Benzyl halides',
	'description': 'R1, R2 = any atom/ group; X = F, Cl, Br, I; A = any aromatic atom',
	'SMARTS': 'a[CX4][F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Benzyl fluorides',
	'description': 'R1, R2 = any atom/ group; A = any aromatic atom',
	'SMARTS': 'a[CX4]F',
	'redux': True,
},
{	
	'name': 'Benzyl chlorides',
	'description': 'R1, R2 = any atom/ group; A = any aromatic atom',
	'SMARTS': 'a[CX4]Cl',
	'redux': True,
},
{	
	'name': 'Benzyl bromides',
	'description': 'R1, R2 = any atom/ group; A = any aromatic atom',
	'SMARTS': 'a[CX4]Br',
	'redux': True,
},
{	
	'name': 'Benzyl iodides',
	'description': 'R1, R2 = any atom/ group; A = any aromatic atom',
	'SMARTS': 'a[CX4]I',
	'redux': True,
},
{	
	'name': 'Alkali metals',
	'description': '',
	'SMARTS': '[#3,#11,#19,#37,#55,#87]',
	'redux': True,
},
{	
	'name': 'Alkaline earth metals',
	'description': '',
	'SMARTS': '[#4,#12,#20,#38,#56,#88]',
	'redux': True,
},
{	
	'name': 'Transition metals',
	'description': '',
	'SMARTS': '[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#72,#73,#74,#75,#76,#77,#78,#79,#80,#104,#105,#106,#107,#108,#109,#110,#111,#112]',
	'redux': True,
},
{	
	'name': 'Lanthanoids',
	'description': '',
	'SMARTS': '[#57,#58,#59,#60,#61,#62,#63,#64,#65,#66,#67,#68,#69,#70,#71]',
	'redux': True,
},
{	
	'name': 'Actinoids',
	'description': '',
	'SMARTS': '[#89,#90,#91,#92,#93,#94,#95,#96,#97,#98,#99,#100,#101,#102,#103]',
	'redux': True,
},
{	
	'name': 'Post-transition metals',
	'description': '',
	'SMARTS': '[#13,#31,#49,#50,#81,#82,#83,#113,#114,#115,#116]',
	'redux': True,
},
{	
	'name': 'Metalloids',
	'description': '',
	'SMARTS': '[#5,#14,#32,#33,#51,#52,#84]',
	'redux': True,
},
{	
	'name': 'Nonmetals',
	'description': '',
	'SMARTS': '[#1,#6,#7,#8,#15,#16,#34,#9,#17,#35,#53,#85,#117,#2,#10,#18,#36,#54,#86,#118]',
	'redux': True,
},
{	
	'name': 'Halogens',
	'description': '',
	'SMARTS': '[#9,#17,#35,#53,#85,#117]',
	'redux': True,
},
{	
	'name': 'Noble gases',
	'description': '',
	'SMARTS': '[#2,#10,#18,#36,#54,#86,#118]',
	'redux': True,
},
{	
	'name': 'Chalcogens (oxygen group)',
	'description': '',
	'SMARTS': '[#8,#16,#34,#52,#84,#116]',
	'redux': True,
},
{	
	'name': 'Pnictogens (nitrogen group)',
	'description': '',
	'SMARTS': '[#7,#15,#33,#51,#83]',
	'redux': True,
},
{	
	'name': 'Tetragens (carbon group)',
	'description': '',
	'SMARTS': '[#6,#14,#32,#50,#8,#114]',
	'redux': True,
},
{	
	'name': 'Thioenolethers',
	'description': 'R1, R2, R3 = H, acyl, alkyl, aryl; R4 = alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[CX3]([SX2][#6&!$([CX3]([OX2])=[OX1,SX1,NX2,C])])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Thioenols',
	'description': 'R1, R2, R3 = H, acyl, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([#1,#6])=[CX3]([SH1])[#1,#6]',
	'redux': True,
},
{	
	'name': 'Enethiodiols',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX3]([#1,#6])([SH1])=[CX3]([#1,#6])([SH1])',
	'redux': True,
},
{	
	'name': 'Alkynyl alcohols',
	'description': 'R = any atom/ group',
	'SMARTS': '[CX2]#[CX2][OH1]',
	'redux': True,
},
{	
	'name': 'Alkynyl thiols',
	'description': 'R = any atom/ group',
	'SMARTS': '[CX2]#[CX2][SH1]',
	'redux': True,
},
{	
	'name': 'Quinones',
	'description': '',
	'SMARTS': '[$([#6X3]1=,:[#6X3]-,:[#6X3](=[OX1])-,:[#6X3]=,:[#6X3]-,:[#6X3]1(=[OX1])),$([#6X3]1(=[OX1])-,:[#6X3](=[OX1])-,:[#6X3]=,:[#6X3]-,:[#6X3]=,:[#6X3]1)]',
	'redux': True,
},
{	
	'name': 'Thioquinones',
	'description': '',
	'SMARTS': '[$([#6X3]1=,:[#6X3]-,:[#6X3](=[SX1])-,:[#6X3]=,:[#6X3]-,:[#6X3]1(=[OX1])),$([#6X3]1(=[OX1])-,:[#6X3](=[OX1,SX1])-,:[#6X3]=,:[#6X3]-,:[#6X3]=,:[#6X3]1)]',
	'redux': True,
},
{	
	'name': 'Three-membered heterocycles (LS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[!#6!#1]1~*~*~1',
	'redux': False,
},
{	
	'name': 'Three-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[!#6!#1]1~[#6]~[#6]~1',
	'redux': False,
},
{	
	'name': 'Saturated three-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon',
	'SMARTS': '[!#6!#1]1-[#6]-[#6]-1',
	'redux': False,
},
{	
	'name': 'Unsaturated three-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[$([!#1!#6]1~[#6]=[#6]~1),$([!#1!#6]1=[#6]~[#6]~1)]',
	'redux': False,
},
{	
	'name': 'Three-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[#6]~1',
	'redux': False,
},
{	
	'name': 'Saturated three-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[#6]-1',
	'redux': False,
},
{	
	'name': 'Unsaturated three-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[$([!#1!#6]1=[!#1!#6]~[#6]~1),$([!#1!#6]1~[!#6!#1]=[#6]~1)]',
	'redux': False,
},

{
	'name': 'Three-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~1',
	'redux': False,
},
{	
	'name': 'Saturated three-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-1',
	'redux': False,
},
{	
	'name': 'Unsaturated three-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[!#1!#6]1~[!#1!#6]=[!#1!#6]~1',
	'redux': False,
},

{
	'name': 'Three-membered heterocycles (HS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[*R1]~[*R1]~1',
	'redux': False,
},

{
	'name': 'Three-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[#6&R1]~[#6&R1]~1',
	'redux': False,
},

{
	'name': 'Saturated three-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[#6&R1]-[#6&R1]-1',
	'redux': False,
},

{
	'name': 'Aziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Oxiranes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Thiiranes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Unsaturated three-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#1!#6;R1]1~[#6R1]=[#6R1]~1),$([!#1!#6;R1]1=[#6R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Azirines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([#7R1]1-[#6R1]=[#6R1]-1),$([#7R1]1=[#6R1]-[#6R1]-1)]',
	'redux': True,
},

{
	'name': '1H-Azirines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#6R1]=[#6R1]-1',
	'redux': True,
},

{
	'name': '2H-Azirines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Oxirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#6R1]=[#6R1]-1',
	'redux': True,
},

{
	'name': 'Thiirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#6R1]=[#6R1]-1',
	'redux': True,
},

{
	'name': 'Three-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated three-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[#6R1]-1',
	'redux': False,
},

{
	'name': 'Diaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#7R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Dioxiranes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#8R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Oxaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#7R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Thiaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#7R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Unsaturated three-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#1!#6;R1]1=[!#1!#6;R1]~[#6R1]~1),$([!#1!#6;R1]1~[!#6!#1;R1]=[#6R1]~1)]',
	'redux': False,
},

{
	'name': '1H-Diazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#7R1]=[#6R1]-1',
	'redux': True,
},

{
	'name': '3H-Diazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#7R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Oxazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#6R1]-[#8R1]-1',
	'redux': True,
},

{
	'name': 'Thiazerenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#6R1]-[#16X2R1]-1',
	'redux': True,
},

{
	'name': 'Three-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~1',
	'redux': False,
},

{
	'name': 'Saturated three-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[!#6!#1;R1]-1',
	'redux': False,
},

{
	'name': 'Triaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#7R1]-[#7R1]-1',
	'redux': True,
},

{
	'name': 'Trioxiranes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#8R1]-[#8R1]-1',
	'redux': True,
},

{
	'name': 'Trithiiranes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#16X2R1]-[#16X2R1]1',
	'redux': True,
},

{
	'name': 'Oxadiaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#7R1]-[#7R1]1',
	'redux': True,
},

{
	'name': 'Dioxaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#8R1]-[#7R1]1',
	'redux': True,
},

{
	'name': 'Thiodiaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#7R1]-[#7R1]1',
	'redux': True,
},

{
	'name': 'Dithiaziridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#16X2R1]-[#7R1]1',
	'redux': True,
},

{
	'name': 'Unsaturated three-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#1!#6;R1]1-[!#6!#1;R1]=[!#6!#1;R1]1',
	'redux': False,
},

{
	'name': '1H-Triazirene (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#7R1]=[#7R1]1',
	'redux': True,
},

{
	'name': 'Oxadiazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#7R1]-[#8R1]1',
	'redux': True,
},

{
	'name': 'Thiodiazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#7R1]-[#16X2R1]1',
	'redux': True,
},

{
	'name': 'Four-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1-[!#6!#1]-[#6]-[#6]-1),$([!#6!#1]1-[#6]-[!#6!#1]-[#6]-1)]',
	'redux': False,
},

{
	'name': 'Unsaturated four-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]=[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]=[#6]~1),$([!#6!#1]1=[#6]~[!#6!#1]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Four-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[#6]-1',
	'redux': False,
},

{
	'name': 'Unsaturated four-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]=[#6]~1)]',
	'redux': False,
},

{
	'name': 'Four-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~1',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[!#6!#1]-1',
	'redux': False,
},

{
	'name': 'Unsaturated four-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1=[!#6!#1]~[!#6!#1]~[!#6!#1]~1',
	'redux': False,
},
{	
	'name': 'Four-membered heterocycles (HS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds',
	'SMARTS': '[!#6!#1;R1]1~[*R1]~[*R1]~[*R1]~1',
	'redux': False,
},

{
	'name': 'Four-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[#6R1]~[#6R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[#6R1]-[#6R1]-[#6R1]-1',
	'redux': False,
},

{
	'name': 'Azetidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#6R1]-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Oxetanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-[#6R1]-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Thietanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1-[#6R1]-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Unsaturated four-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]=[#6R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Azetines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([#7R1]1=[#6R1]-[#6R1]-[#6R1]-1),$([#7R1]1-[#6R1]=[#6R1]-[#6R1]-1)]',
	'redux': True,
},

{
	'name': '1-Azetines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#6R1]-[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': '2-Azetines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1-[#6R1]=[#6R1]-[#6R1]-1',
	'redux': True,
},

{
	'name': 'Azetes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#7R1]1=[#6R1][#6R1]=[#6R1]1',
	'redux': True,
},

{
	'name': 'Oxetenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1[#6R1]=[#6R1][#6R1]1',
	'redux': True,
},

{
	'name': 'Thietenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16X2R1]1[#6R1]=[#6R1][#6R1]1',
	'redux': True,
},

{
	'name': 'Four-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1;R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1-[!#6!#1;R1]-[#6R1]-[#6R1]-1),$([!#6!#1;R1]1-[#6R1]-[!#6!#1;R1]-[#6R1]-1)]',
	'redux': False,
},

{
	'name': 'Diazetidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[CR1][CR1][NR1]1),$([NR1]1[CR1][NR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Diazetidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][CR1][NR1]1',
	'redux': True,
},

{
	'name': '1,3-Diazetidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][NR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dioxetanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][CR1][CR1]1),$([OR1]1[CR1][OR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Dioxetanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dioxetanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][OR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dithietanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SX2R1]1[SX2R1][CR1][CR1]1),$([SX2R1]1[CR1][SX2R1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Dithietanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SX2R1]1[SX2R1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dithietanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SX2R1]1[CR1][SX2R1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated four-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]=[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]=[#6R1]~1),$([!#6!#1;R1]1=[#6R1]~[!#6!#1;R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Diazetines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([#7]1=,:[#7]-,:[#6]-,=,:[#6]-,:1),$([#7]1-,:[#7]=,:[#6]-,:[#6]-,=,:1),$([#7]1-,:[#7]-,:[#6]=,:[#6]-,=,:1),$([#7]1=,:[#6]-,:[#7]-,=,:[#6]-,:1)]',
	'redux': True,
},

{
	'name': 'Dioxetenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#8R1]1-,:[#8R1]-,:[#6R1]=,:[#6R1]-,:1',
	'redux': True,
},

{
	'name': 'Dithietenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[#16R1]1-,:[#16R1]-,:[#6R1]=,:[#6R1]-,:1',
	'redux': True,
},

{
	'name': 'Four-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[#6]~[#6]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated four-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[#6]-[#6]-[#6]-1',
	'redux': False,
},

{
	'name': 'Unsaturated four-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[#6]~[#6]~[#6]~1),$([!#6!#1]1~[#6]=[#6]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Diazirenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([#7R1]1-[#7R1]=[#6R1]-1),$([#7R1]1=[#7R1]-[#6R1]-1)]',
	'redux': True,
},
{	
	'name': 'Unsaturated three-membered heterocycles with two heteroatoms (LS)',
	'description': '',
	'SMARTS': '[$([!#1!#6]1=[!#1!#6]-[#6]1),$([!#1!#6]1-[!#6!#1]=[#6]1)]',
	'redux': False,
},

{
	'name': 'Four-membered heterocycles (LS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~*~*~*~1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles (LS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~*~*~*~*~1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[#6]~[#6]~[#6]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1CCCC1',
	'redux': False,
},

{
	'name': 'Unsaturated five-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=CCCC1),$([!#6!#1]1C=CCC1),$([!#6!#1]1CC=CC1),$([!#6!#1]1=CC=CC1);!$(n1cccc1)]',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!#6]1cccc1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1[!#6!#1]CCC1),$([!#6!#1]1C[!#6!#1]CC1)]',
	'redux': False,
},

{
	'name': 'Unsaturated five-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]C-,=CC1),$([!#6!#1]1[!#6!#1]=CCC1),$([!#6!#1]1[!#6!#1]C=CC1),$([!#6!#1]1=C[!#6!#1]-,=CC1),$([!#6!#1]1C[!#6!#1]CC=1),$([!#6!#1]1C[!#6!#1]C=C1)]',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([a!#6]1[a!#6]ccc1),$([a!#6]1c[a!#6]cc1)]',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1[!#6!#1][!#6!#1]CC1),$([!#6!#1]1[!#6!#1]C[!#6!#1]C1)]',
	'redux': False,
},

{
	'name': 'Unsaturated five-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1][!#6!#1]-,=CC1),$([!#6!#1]1[!#6!#1][!#6!#1]CC=1),$([!#6!#1]1[!#6!#1][!#6!#1]C=C1),$([!#6!#1]1=[!#6!#1]C[!#6!#1]-,=C1),$([!#6!#1]1[!#6!#1]C[!#6!#1]C=1),$([!#6!#1]1[!#6!#1]C[!#6!#1]=C1)]',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([a!#6]1[a!#6][a!#6]cc1),$([a!#6]1[a!#6]c[a!#6]c1)]',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[!#6!#1]C1',
	'redux': False,
},

{
	'name': 'Unsaturated five-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1][!#6!#1]-,=[!#6!#1][#6]1),$([!#6!#1]1[!#6!#1]=[!#6!#1][!#6!#1][#6]1),$([!#6!#1]1[!#6!#1][!#6!#1][!#6!#1][#6]=1)]',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!#6]1[a!#6][a!#6][a!#6]c1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[!#6!#1]-[!#6!#1]-1',
	'redux': False,
},

{
	'name': 'Unsaturated five-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1=[!#6!#1]-[!#6!#1]-[!#6!#1]-[!#6!#1]-1',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!#6]1[a!#6][a!#6][a!#6][a!#6]1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles (HS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[*R1]~[*R1]~[*R1]~[*R1]~1',
	'redux': False,
},

{
	'name': 'Five-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[#6R1]~[#6R1]~[#6R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1[CR1][CR1][CR1][CR1]1',
	'redux': False,
},

{
	'name': 'Pyrrolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetrahydrofurans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetrahydrothiophenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated five-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[CR1][CR1][CR1][CR1]1),$([!#6!#1;R1]1[CR1]=[CR1][CR1][CR1]1),$([!#6!#1;R1]1[CR1][CR1]=[CR1][CR1]1),$([!#6!#1;R1]1=[CR1][CR1]=[CR1][CR1]1);!$(n1cccc1)]',
	'redux': False,
},

{
	'name': 'Pyrrolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[CR1][CR1][CR1][CR1]1),$([NR1]1[CR1]=[CR1][CR1][CR1]1),$([NR1]1[CR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1-Pyrrolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1=[CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '2-Pyrrolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1]=[CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '3-Pyrrolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][CR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dihydrofurans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[CR1]=[CR1][CR1][CR1]1),$([OR1]1[CR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '2,3-Dihydrofurans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1]=[CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '2,5-Dihydrofurans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][CR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dihydrothiophenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[CR1]=[CR1][CR1][CR1]1),$([SR1]1[CR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '2,3-Dihydrothiophenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1]=[CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '2,5-Dihydrothiophenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1][CR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Aromatic five-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!#6;R1]1[cR1][cR1][cR1][cR1]1',
	'redux': False,
},

{
	'name': 'Pyrroles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[cR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Furans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[cR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Thiophenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[sR1]1[cR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Indoles',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[cR1]1[cR1][cR1][cR1][cR2]2[nR1][cR1][cR1][cR2]12',
	'redux': True,
},

{
	'name': 'Benzofurans',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[cR1]1[cR1][cR1][cR1][cR2]2[oR1][cR1][cR1][cR2]12',
	'redux': True,
},

{
	'name': 'Benzothiophenes',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[cR1]1[cR1][cR1][cR1][cR2]2[sR1][cR1][cR1][cR2]12',
	'redux': True,
},

{
	'name': 'Five-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1[!#6!#1;R1][CR1][CR1][CR1]1),$([!#6!#1;R1]1[CR1][!#6!#1;R1][CR1][CR1]1)]',
	'redux': False,
},

{
	'name': 'Diazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[NR1][CR1][CR1][CR1]1),$([NR1]1[CR1][NR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Pyrazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Imidazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][NR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][CR1][CR1][CR1]1),$([OR1]1[CR1][OR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Dioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][OR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Dithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][CR1][CR1][CR1]1),$([SR1]1[CR1][SR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Dithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1][SR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated five-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1][CR1]-,=[CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1]=[CR1][CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][CR1]=[CR1][CR1]1),$([!#6!#1;R1]1=[CR1][!#6!#1;R1]-,=[CR1][CR1]1),$([!#6!#1;R1]1[CR1][!#6!#1;R1][CR1][CR1]=1),$([!#6!#1;R1]1[CR1][!#6!#1;R1][CR1]=[CR1]1)]',
	'redux': False,
},

{
	'name': 'Pyrazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][CR1][CR1][CR1]1),$([NR1]1[NR1]=[CR1][CR1][CR1]1),$([NR1]1[NR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1-Pyrazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1=[NR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '2-Pyrazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1]=[CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Imidazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[CR1][NR1][CR1][CR1]1),$([NR1]1[CR1][NR1]=[CR1][CR1]1),$([NR1]1[CR1][NR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': '2-Imidazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1=[CR1][NR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '3-Imidazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][NR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': '4-Imidazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][NR1][CR1]=[CR1]1',
	'redux': True,
},

{
	'name': 'Dioxolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][CR1]=[CR1][CR1]1),$([OR1]1[CR1][OR1][CR1]=[CR1]1)]',
	'redux': True,
},
{
	'name': '1,2-Dioxolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][CR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dioxolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][OR1][CR1]=[CR1]1',
	'redux': True,
},

{
	'name': 'Dithiolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][CR1]=[CR1][CR1]1),$([SR1]1[CR1][SR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Dithiolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][CR1]=[CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,3-Dithiolenes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1][SR1][CR1]=[CR1]1',
	'redux': True,
},

{
	'name': 'Aromatic five-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([a!#6;R1]1[a!#6;R1][cR1][cR1][cR1]1),$([a!#6;R1]1[cR1][a!#6;R1][cR1][cR1]1)]',
	'redux': False,
},

{
	'name': 'Diazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([nR1]1[nR1][cR1][cR1][cR1]1),$([nR1]1[cR1][nR1][cR1][cR1]1)]',
	'redux': True,
},

{
	'name': 'Pyrazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Imidazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[cR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Oxazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([oR1]1[nR1][cR1][cR1][cR1]1),$([oR1]1[cR1][nR1][cR1][cR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Oxazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': '1,3-Oxazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[cR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Thiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([sR1]1[nR1][cR1][cR1][cR1]1),$([sR1]1[cR1][nR1][cR1][cR1]1)]',
	'redux': True,
},

{
	'name': '1,2-Thiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[sR1]1[nR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': '1,3-Thiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[sR1]1[cR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Five-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][CR1][!#6!#1;R1][CR1]1)]',
	'redux': False,
},

{
	'name': 'Triazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[NR1][NR1][CR1][CR1]1),$([NR1]1[NR1][CR1][NR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Triazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][NR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,2,4-Triazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][CR1][NR1][CR1]1',
	'redux': True,
},

{
	'name': 'Trioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][OR1][CR1][CR1]1),$([OR1]1[OR1][CR1][OR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Trioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][OR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,2,4-Trioxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][CR1][OR1][CR1]1',
	'redux': True,
},

{
	'name': 'Trithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][SR1][CR1][CR1]1),$([SR1]1[SR1][CR1][SR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Trithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][SR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': '1,2,4-Trithiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][CR1][SR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated five-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1][!#6!#1;R1]-,=[CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][CR1][CR1]=1),$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][CR1]=[CR1]1),$([!#6!#1;R1]1=[!#6!#1;R1][CR1][!#6!#1;R1]-,=[CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][CR1][!#6!#1;R1][CR1]=1),$([!#6!#1;R1]1[!#6!#1;R1][CR1][!#6!#1;R1]=[CR1]1)]',
	'redux': False,
},

{
	'name': '1,2,3-Triazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][NR1][CR1][CR1]1),$([NR1]1[NR1][NR1][CR1][CR1]=1),$([NR1]1[NR1][NR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': '1,2,4-Triazolines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][CR1][NR1][CR1]1),$([NR1]1[NR1]=[CR1][NR1][CR1]1),$([NR1]1[NR1][CR1]=[NR1][CR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Trioxoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][OR1][CR1]=[CR1]1',
	'redux': True,
},

{
	'name': '1,2,3-Trithioles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][SR1][CR1]=[CR1]1',
	'redux': True,
},

{
	'name': 'Aromatic five-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([a!#6;R1]1[a!#6;R1][a!#6;R1][cR1][cR1]1),$([a!#6;R1]1[a!#6;R1][cR1][a!#6;R1][cR1]1)]',
	'redux': False,
},

{
	'name': 'Triazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([nR1]1[nR1][nR1][cR1][cR1]1),$([nR1]1[nR1][cR1][nR1][cR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Triazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': '1,2,4-Triazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][cR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': 'Oxadiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([oR1]1[nR1][nR1][cR1][cR1]1),$([oR1]1[nR1][cR1][nR1][cR1]1),$([oR1]1[nR1][cR1][cR1][nR1]1),$([oR1]1[cR1][nR1][nR1][cR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3-Oxadiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': '1,2,4-Oxadiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][cR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': '1,2,5-Oxadiazoles (furazanes) (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][cR1][cR1][nR1]1',
	'redux': True,
},

{
	'name': '1,3,4-Oxadiazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[cR1][nR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': 'Five-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1][CR1]1',
	'redux': False,
},

{
	'name': 'Tetrazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][NR1][NR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetraoxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][OR1][OR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetrathiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][SR1][SR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated five-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1][!#6!#1;R1]-,=[!#6!#1;R1][#6R1]1),$([!#6!#1;R1]1[!#6!#1;R1]=[!#6!#1;R1][!#6!#1;R1][#6R1]1),$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][!#6!#1;R1][#6R1]=1)]',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!#6;R1]1[a!#6;R1][a!#6;R1][a!#6;R1][cR1]1',
	'redux': False,
},

{
	'name': 'Tetrazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][nR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': 'Oxatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([oR1]1[nR1][nR1][nR1][cR1]1),$([oR1]1[nR1][nR1][cR1][nR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3,4-Oxatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][nR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': '1,2,3,5-Oxatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[oR1]1[nR1][nR1][cR1][nR1]1',
	'redux': True,
},

{
	'name': 'Thiatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([sR1]1[nR1][nR1][nR1][cR1]1),$([sR1]1[nR1][nR1][cR1][nR1]1)]',
	'redux': True,
},

{
	'name': '1,2,3,4-Thiatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[sR1]1[nR1][nR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': '1,2,3,5-Thiatriazoles (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[sR1]1[nR1][nR1][cR1][nR1]1',
	'redux': True,
},

{
	'name': 'Five-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~1',
	'redux': False,
},

{
	'name': 'Saturated five-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-1',
	'redux': False,
},

{
	'name': 'Pentazolidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][NR1][NR1][NR1]1',
	'redux': True,
},

{
	'name': 'Pentaoxolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][OR1][OR1][OR1]1',
	'redux': True,
},

{
	'name': 'Pentathiolanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][SR1][SR1][SR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated five-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1=[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-1',
	'redux': False,
},

{
	'name': 'Aromatic five-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!#6;R1]1[a!#6;R1][a!#6;R1][a!#6;R1][a!#6;R1]1',
	'redux': False,
},

{
	'name': 'Pentazoles',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][nR1][nR1][nR1]1',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles (LS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~*~*~*~*~*~1',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[#6]~[#6]~[#6]~[#6]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1CCCCC1',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#1!#6]1=C~C~C~C~C~1),$([!#1!#6]1~C=C~C~C~C~1),$([!#1!#6]1~C~C=C~C~C~1);!$([a!c]1ccccc1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with one heteroatom (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!c]1ccccc1',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[#6]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[#6]~[!#6!#1]~[#6]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1[!#6!#1]CCCC1),$([!#6!#1]1C[!#6!#1]CCC1),$([!#6!#1]1CC[!#6!#1]CC1)]',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[#6]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]=[#6]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]=[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[#6]=[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~[#6]~[#6]=1),$([!#6!#1]1=[#6]~[!#6!#1]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~[#6]=[#6]~1),$([!#6!#1]1=[#6]~[#6]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[#6]=[#6]~[!#6!#1]~[#6]~[#6]~1);!$([a!c]1[a!c]cccc1);!$([a!c]1c[a!c]ccc1);!$([a!c]1cc[a!c]cc1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with two heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([a!c]1[a!c]cccc1),$([a!c]1c[a!c]ccc1),$([a!c]1cc[a!c]cc1)]',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[#6]~[!#6!#1]~[#6]~[!#6!#1]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1[!#6!#1][!#6!#1]CCC1),$([!#6!#1]1[!#6!#1]C[!#6!#1]CC1),$([!#6!#1]1C[!#6!#1]C[!#6!#1]C1)]',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]=[#6]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]=[#6]~[#6]~1),$([!#6!#1]1=[!#6!#1]~[#6]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]=[#6]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]=[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]=[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]~[#6]=[#6]~1),$([!#6!#1]1=[#6]~[!#6!#1]~[#6]~[!#6!#1]~[#6]~1);!$([a!c]1[a!c][a!c]ccc1);!$([a!c]1[a!c]c[a!c]cc1);!$([a!c]1c[a!c]c[a!c]c1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with three heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([a!c]1[a!c][a!c]ccc1),$([a!c]1[a!c]c[a!c]cc1),$([a!c]1c[a!c]c[a!c]c1)]',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~[!#6!#1]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]~[!#6!#1]~[#6]~1)]',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1[!#6!#1][!#6!#1][!#6!#1]CC1),$([!#6!#1]1[!#6!#1][!#6!#1]C[!#6!#1]C1),$([!#6!#1]1[!#6!#1]C[!#6!#1][!#6!#1]C1)]',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]=[!#6!#1]~[!#6!#1]~[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]=[#6]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]=[#6]~1),$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[#6]~[!#6!#1]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~[!#6!#1]~[#6]=1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[#6]~[!#6!#1]=[#6]~1),$([!#6!#1]1~[!#6!#1]~[#6]~[!#6!#1]~[!#6!#1]~[#6]=1),$([!#6!#1]1=[!#6!#1]~[#6]~[!#6!#1]~[!#6!#1]~[#6]~1);!$([a!c]1[a!c][a!c][a!c]cc1);!$([a!c]1[a!c][a!c]c[a!c]c1);!$([a!c]1[a!c]c[a!c][a!c]c1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with four heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([a!c]1[a!c][a!c][a!c]cc1),$([a!c]1[a!c][a!c]c[a!c]c1),$([a!c]1[a!c]c[a!c][a!c]c1)]',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[!#6!#1]-[!#6!#1]-C-1',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~1),$([!#6!#1]1~[!#6!#1]=[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]~1),$([!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~[#6]=1);!$([a!c]1[a!c][a!c][a!c][a!c]c1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with five heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!c]1[a!c][a!c][a!c][a!c]c1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with six heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1-[!#6!#1]-[!#6!#1]-[!#6!#1]-[!#6!#1]-[!#6!#1]-1',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with six heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[$([!#6!#1]1=[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~1);!$([a!c]1[a!c][a!c][a!c][a!c][a!c]1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with six heteroatoms (LS)',
	'description': 'A = any aromatic atom except carbon; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[a!c]1[a!c][a!c][a!c][a!c][a!c]1',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles (HS)',
	'description': 'A = any atom except carbon; A1 = any atom; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[*R1]~[*R1]~[*R1]~[*R1]~[*R1]~1',
	'redux': False,
},

{
	'name': 'Six-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[#6R1]~[#6R1]~[#6R1]~[#6R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1[CR1][CR1][CR1][CR1][CR1]1',
	'redux': False,
},

{
	'name': 'Piperidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetrahydropyrans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[CR1][CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Tetrahydrothiopyrans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1][CR1][CR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated six-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#1!#6;R1]1=[CR1]~[CR1]~[CR1]~[CR1]~[CR1]~1),$([!#1!#6;R1]1~[CR1]=[CR1]~[CR1]~[CR1]~[CR1]~1),$([!#1!#6;R1]1~[CR1]~[CR1]=[CR1]~[CR1]~[CR1]~1);!$([a!c]1ccccc1)]',
	'redux': False,
},

{
	'name': 'Tetrahydropyridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[CR1][CR1][CR1][CR1][CR1]1),$([NR1]1[CR1]=[CR1][CR1][CR1][CR1]1),$([NR1]1[CR1][CR1]=[CR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Dihydropyridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[CR1][CR1]-,=[CR1]-,=[CR1][CR1]1),$([NR1]1[CR1]=[CR1][CR1][CR1]-,=[CR1]-,=1),$([NR1]1[CR1]=[CR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Pyrans (HS)',
	'description': '; Low specifisity (HS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[OR1]1[CR1]=[CR1][CR1]-,=[CR1]-,=[CR1]1',
	'redux': True,
},

{
	'name': 'Thiopyrans (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[CR1]=[CR1][CR1]-,=[CR1]-,=[CR1]1',
	'redux': True,
},

{
	'name': 'Aromatic six-membered heterocycles with one heteroatom (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!c;R1]1[cR1][cR1][cR1][cR1][cR1]1',
	'redux': False,
},

{
	'name': 'Pyridines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[cR1][cR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1]~[#6R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[#6R1]~[!#6!#1]~[#6R1]~[#6R1]~1)]',
	'redux': False,
},
{
	'name': 'Saturated six-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1[!#6!#1;R1][CR1][CR1][CR1][CR1]1),$([!#6!#1;R1]1[CR1][!#6!#1;R1][CR1][CR1][CR1]1),$([!#6!#1;R1]1[CR1][CR1][!#6!#1;R1][CR1][CR1]1)]',
	'redux': False,
},

{
	'name': 'Hexahydrodiazines',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[NR1][CR1][CR1][CR1][CR1]1),$([NR1]1[CR1][NR1][CR1][CR1][CR1]1),$([NR1]1[CR1][CR1][NR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Dioxanes',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][CR1][CR1][CR1][CR1]1),$([OR1]1[CR1][OR1][CR1][CR1][CR1]1),$([OR1]1[CR1][CR1][OR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Dithianes',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][CR1][CR1][CR1][CR1]1),$([SR1]1[CR1][SR1][CR1][CR1][CR1]1),$([SR1]1[CR1][CR1][SR1][CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Morpholines',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[CR1][CR1][OR1][CR1][CR1]1',
	'redux': True,
},

{
	'name': 'Unsaturated six-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]=[#6R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]=[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[#6R1]=[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]=1),$([!#6!#1;R1]1=[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]=[#6R1]~1),$([!#6!#1;R1]1=[#6R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]=[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1);!$([a!c]1[a!c]cccc1);!$([a!c]1c[a!c]ccc1);!$([a!c]1cc[a!c]cc1)]',
	'redux': False,
},

{
	'name': 'Tetrahydrodiazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][CR1][CR1][CR1][CR1]1),$([NR1]1[NR1]=[CR1][CR1][CR1][CR1]1),$([NR1]1[NR1][CR1]=[CR1][CR1][CR1]1),$([NR1]1[NR1][CR1][CR1]=[CR1][CR1]1),$([NR1]1=[CR1][CR1][NR1][CR1][CR1]1),$([NR1]1[CR1]=[CR1][NR1][CR1][CR1]1),$([NR1]1=[CR1][NR1][CR1][CR1][CR1]1),$([NR1]1[CR1][NR1]=[CR1][CR1][CR1]1),$([NR1]1[CR1][NR1][CR1]=[CR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Dihydrodiazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][CR1]=[CR1][CR1][CR1]1),$([NR1]1=[NR1][CR1]=[CR1][CR1][CR1]1),$([NR1]1[NR1][CR1]=[CR1][CR1]=[CR1]1),$([NR1]1[NR1]=[CR1][CR1]=[CR1][CR1]1),$([NR1]1[NR1]=[CR1][CR1][CR1]=[CR1]1),$([NR1]1[NR1]=[CR1][CR1][CR1][CR1]=1),$([NR1]1[CR1][NR1]=[CR1][CR1][CR1]=1),$([NR1]1=[CR1][NR1]=[CR1][CR1][CR1]1),$([NR1]1=[CR1][NR1][CR1]=[CR1][CR1]1),$([NR1]1=[CR1][NR1][CR1][CR1]=[CR1]1),$([NR1]1[CR1][NR1]=[CR1][CR1]=[CR1]1),$([NR1]1=[CR1][CR1]=[NR1][CR1][CR1]1),$([NR1]1=[CR1][CR1][NR1]=[CR1][CR1]1),$([NR1]1=[CR1][CR1][NR1][CR1]=[CR1]1),$([NR1]1[CR1]=[CR1][NR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': 'Dioxines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][CR1]=[CR1][CR1]=[CR1]1),$([OR1]1[CR1]=[CR1][OR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': 'Aromatic six-membered heterocycles with two heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([a!c;R1]1[a!c;R1][cR1][cR1][cR1][cR1]1),$([a!c;R1]1[cR1][a!c;R1][cR1][cR1][cR1]1),$([a!c;R1]1[cR1][cR1][a!c;R1][cR1][cR1]1)]',
	'redux': False,
},

{
	'name': 'Pyridazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][cR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Pyrimidines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[cR1][nR1][cR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Pyrazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[cR1][cR1][nR1][cR1][cR1]1',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[#6R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][CR1][CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][CR1][!#6!#1;R1][CR1][CR1]1),$([!#6!#1;R1]1[CR1][!#6!#1;R1][CR1][!#6!#1;R1][CR1]1)]',
	'redux': False,
},

{
	'name': 'Hexahydrotriazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[NR1][NR1][CR1][CR1][CR1]1),$([NR1]1[NR1][CR1][NR1][CR1][CR1]1),$([NR1]1[CR1][NR1][CR1][NR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Trioxanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][OR1][CR1][CR1][CR1]1),$([OR1]1[OR1][CR1][OR1][CR1][CR1]1),$([OR1]1[CR1][OR1][CR1][OR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Trithianes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][SR1][CR1][CR1][CR1]1),$([SR1]1[SR1][CR1][SR1][CR1][CR1]1),$([SR1]1[CR1][SR1][CR1][SR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Unsaturated six-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]=[#6R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]=[#6R1]~[#6R1]~1),$([!#6!#1;R1]1=[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]=[#6R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]=[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]=[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]=[#6R1]~1),$([!#6!#1;R1]1=[#6R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~1);!$([a!c]1[a!c][a!c]ccc1);!$([a!c]1[a!c]c[a!c]cc1);!$([a!c]1c[a!c]c[a!c]c1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with three heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([a!c;R1]1[a!c;R1][a!c;R1][cR1][cR1][cR1]1),$([a!c;R1]1[a!c;R1][cR1][a!c;R1][cR1][cR1]1),$([a!c;R1]1[cR1][a!c;R1][cR1][a!c;R1][cR1]1)]',
	'redux': False,
},

{
	'name': 'Triazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([nR1]1[nR1][nR1][cR1][cR1][cR1]1),$([nR1]1[nR1][cR1][nR1][cR1][cR1]1),$([nR1]1[cR1][nR1][cR1][nR1][cR1]1)]',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1)]',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][!#6!#1;R1][CR1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][!#6!#1;R1][CR1][!#6!#1;R1][CR1]1),$([!#6!#1;R1]1[!#6!#1;R1][CR1][!#6!#1;R1][!#6!#1;R1][CR1]1)]',
	'redux': False,
},

{
	'name': 'Hexahydrotetrazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1[NR1][NR1][NR1][CR1][CR1]1),$([NR1]1[NR1][CR1][NR1][NR1][CR1]1),$([NR1]1[NR1][NR1][CR1][NR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Tetroxanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([OR1]1[OR1][OR1][OR1][CR1][CR1]1),$([OR1]1[OR1][CR1][OR1][OR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Tetrathianes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][SR1][SR1][CR1][CR1]1),$([SR1]1[SR1][CR1][SR1][SR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Unsaturated six-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]=[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]=[#6R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]=[#6R1]~1),$([!#6!#1;R1]1=[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[#6R1]=1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]=[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]=1),$([!#6!#1;R1]1=[!#6!#1;R1]~[#6R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1);!$([a!c]1[a!c][a!c][a!c]cc1);!$([a!c]1[a!c][a!c]c[a!c]c1);!$([a!c]1[a!c]c[a!c][a!c]c1)]',
	'redux': False,
},

{
	'name': 'Tetrahydrotetrazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][NR1][NR1][CR1][CR1]1),$([NR1]1[NR1]=[NR1][NR1][CR1][CR1]1),$([NR1]1[NR1][NR1][NR1][CR1][CR1]=1),$([NR1]1[NR1][NR1][NR1][CR1]=[CR1]1),$([NR1]1[NR1][NR1][CR1][NR1][CR1]=1),$([NR1]1[NR1][NR1][CR1][NR1]=[CR1]1),$([NR1]1[NR1]=[NR1][CR1][NR1][CR1]1),$([NR1]1[NR1]=[CR1][NR1][NR1][CR1]1),$([NR1]1=[NR1][CR1][NR1][NR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Aromatic six-membered heterocycles with four heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([a!c;R1]1[a!c;R1][a!c;R1][a!c;R1][cR1][cR1]1),$([a!c;R1]1[a!c;R1][a!c;R1][cR1][a!c;R1][cR1]1),$([a!c;R1]1[a!c;R1][cR1][a!c;R1][a!c;R1][cR1]1)]',
	'redux': False,
},

{
	'name': 'Tetrazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([nR1]1[nR1][nR1][nR1][cR1][cR1]1),$([nR1]1[nR1][nR1][cR1][nR1][cR1]1),$([nR1]1[nR1][cR1][nR1][nR1][cR1]1)]',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[CR1]-1',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]=[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]~1),$([!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[#6R1]=1);!$([a!c]1[a!c][a!c][a!c][a!c]c1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with five heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!c;R1]1[a!c;R1][a!c;R1][a!c;R1][a!c;R1][cR1]1',
	'redux': False,
},

{
	'name': 'Pentazines',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][nR1][nR1][nR1][cR1]1',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with six heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~1',
	'redux': False,
},

{
	'name': 'Saturated six-membered heterocycles with six heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[!#6!#1;R1]1-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-[!#6!#1;R1]-1',
	'redux': False,
},

{
	'name': 'Unsaturated six-membered heterocycles with six heteroatoms (HS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([!#6!#1;R1]1=[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~[!#6!#1;R1]~1);!$([a!c]1[a!c][a!c][a!c][a!c][a!c]1)]',
	'redux': False,
},

{
	'name': 'Aromatic six-membered heterocycles with six heteroatoms (HS)',
	'description': 'A = any aromatic atom except carbon; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[a!c;R1]1[a!c;R1][a!c;R1][a!c;R1][a!c;R1][a!c;R1]1',
	'redux': False,
},

{
	'name': 'Hexazines',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[nR1]1[nR1][nR1][nR1][nR1][nR1]1',
	'redux': True,
},

{
	'name': 'Six-membered heterocycles with six heteroatoms (LS)',
	'description': 'A = any atom except carbon; a dashed line indicates any type of covalent bonds; Low specifisity (LS) pattern matches any chemicals that include depicted heterocyclic moiety (fusion with other ring(s) are allowed)',
	'SMARTS': '[!#6!#1]1~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~[!#6!#1]~1',
	'redux': False,
},

{
	'name': 'Dithiines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([SR1]1[SR1][CR1]=[CR1][CR1]=[CR1]1),$([SR1]1[CR1]=[CR1][SR1][CR1]=[CR1]1)]',
	'redux': True,
},

{
	'name': 'Dihydrotetrazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[$([NR1]1=[NR1][NR1][NR1][CR1]=[CR1]1),$([NR1]1[NR1]=[NR1][NR1][CR1]=[CR1]1),$([NR1]1[NR1]=[NR1][NR1]=[CR1][CR1]1),$([NR1]1[NR1][NR1]=[NR1][CR1][CR1]=1),$([NR1]1[NR1][NR1][NR1]=[CR1][CR1]=1),$([NR1]1=[NR1][NR1]=[NR1][CR1][CR1]1),$([NR1]1=[NR1][NR1][CR1][NR1]=[CR1]1),$([NR1]1[NR1]=[NR1][CR1][NR1]=[CR1]1),$([NR1]1[NR1][NR1]=[CR1][NR1]=[CR1]1),$([NR1]1[NR1]=[NR1][CR1][NR1][CR1]=1),$([NR1]1[NR1][NR1]=[CR1][NR1][CR1]=1),$([NR1]1=[NR1][NR1]=[CR1][NR1][CR1]1),$([NR1]1=[NR1][NR1][CR1]=[NR1][CR1]1),$([NR1]1[NR1]=[CR1][NR1][NR1][CR1]=1),$([NR1]1[NR1][CR1]=[NR1][NR1][CR1]=1),$([NR1]1[NR1][CR1][NR1]=[NR1][CR1]=1),$([NR1]1=[NR1][CR1]=[NR1][NR1][CR1]1),$([NR1]1=[NR1][CR1][NR1]=[NR1][CR1]1)]',
	'redux': True,
},

{
	'name': 'Hexahydropentazines (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[NR1]1[NR1][NR1][NR1][NR1][CR1]1',
	'redux': True,
},

{
	'name': 'Pentoxanes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[OR1]1[OR1][OR1][OR1][OR1][CR1]1',
	'redux': True,
},

{
	'name': 'Pentathianes (HS)',
	'description': '; High specificity (HS) pattern matches chemicals that include exact heterocyclic moiety as in the depiction (fusion with other ring(s) are not allowed)',
	'SMARTS': '[SR1]1[SR1][SR1][SR1][SR1][CR1]1',
	'redux': True,
},
{	
	'name': 'Alkanes',
	'description': 'Compounds consist only of hydrogen and carbon atoms, all bonds are single bonds, carbon atoms can be joined in a chain or cyclic structures',
	'SMARTS': 'NOT [!#6!#1] AND NOT [$([#6]=[#6]),$([#6]#[#6]),$(cc)]',
	'redux': True,
},
{	
	'name': 'Alkenes',
	'description': 'R1, R2, R3, R4 = H, alkyl, aryl',
	'SMARTS': '[!#8!#7!#16][CX3]([!#8!#7!#16])=[CX3]([!#8!#7!#16])[!#8!#7!#16] AND NOT [$([CX3]=[CX2]=[CX3]),$([CX3]=[CX3][CX3]=[CX3]),$([CX3]=[CX3][CX4][CX3]=[CX3])]',
	'redux': True,
},
{	
	'name': 'Cumulated alkadienes (1,2-alkadienes)',
	'description': 'R1, R2, R3, R4 = any atom/ group',
	'SMARTS': '[CX3]=[CX2]=[CX3]',
	'redux': True,
},
{	
	'name': 'Allyl halides',
	'description': 'R1, R2, R3, R4, R5 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4][F,Cl,Br,I]',
	'redux': True,
},
{	
	'name': 'Allyl fluorides',
	'description': 'R1, R2, R3, R4, R5 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4]F',
	'redux': True,
},
{	
	'name': 'Allyl chlorides',
	'description': 'R1, R2, R3, R4, R5 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4]Cl',
	'redux': True,
},
{	
	'name': 'Allyl bromides',
	'description': 'R1, R2, R3, R4, R5 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4]Br',
	'redux': True,
},
{	
	'name': 'Allyl iodides',
	'description': 'R1, R2, R3, R4, R5 = any atom/ group',
	'SMARTS': '[CX3]=[CX3][CX4]I',
	'redux': True,
},
{	
	'name': 'Vinyl fluorides',
	'description': 'R1, R2, R3 = any atom/ group',
	'SMARTS': 'F[CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Vinyl chlorides',
	'description': 'R1, R2, R3 = any atom/ group',
	'SMARTS': 'Cl[CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Vinyl bromides',
	'description': 'R1, R2, R3 = any atom/ group',
	'SMARTS': 'Br[CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Vinyl iodides',
	'description': 'R1, R2, R3 = any atom/ group',
	'SMARTS': 'I[CX3]=[CX3]',
	'redux': True,
},
{	
	'name': 'Alkynyl fluorides',
	'description': 'R = any atom/ group',
	'SMARTS': 'F[CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Alkynyl chlorides',
	'description': 'R = any atom/ group',
	'SMARTS': 'Cl[CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Alkynyl iodides',
	'description': 'R = any atom/ group',
	'SMARTS': 'I[CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Alkynyl bromides',
	'description': 'R = any atom/ group',
	'SMARTS': 'Br[CX2]#[CX2]',
	'redux': True,
},
{	
	'name': 'Dialkylthioethers',
	'description': 'R1, R2 = alkyl',
	'SMARTS': '[CX4&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))][SX2][CX4&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))]',
	'redux': True,
},
{	
	'name': 'Alkylarylthioethers',
	'description': 'R1 = alkyl, R2 = aryl',
	'SMARTS': '[CX4&!$([CX4]([SX2])([#7,O,S,F,Cl,Br,I,P]))][SX2]c',
	'redux': True,
},
{	
	'name': 'Diarylthioethers',
	'description': 'R1, R2 = aryl',
	'SMARTS': 'c[SX2&!$([SX2r3])]c',
	'redux': True,
},
{	
	'name': '1,2 - Diphenols',
	'description': '',
	'SMARTS': 'c([OH1])c([OH1])',
	'redux': True,
},
{	
	'name': '1,3 - Diphenols',
	'description': '',
	'SMARTS': 'c([OH1])[aR1]c([OH1])',
	'redux': True,
},
{	
	'name': '1,4 - Diphenols',
	'description': '',
	'SMARTS': 'c1([OH1])aac([OH1])aa1',
	'redux': True,
},
{	
	'name': '1,2 - Dithiophenols',
	'description': '',
	'SMARTS': 'c([SH1])c([SH1])',
	'redux': True,
},
{	
	'name': '1,3 - Dithiophenols',
	'description': '',
	'SMARTS': 'c([SH1])[aR1]c([SH1])',
	'redux': True,
},
{	
	'name': '1,4 - Dithiophenols',
	'description': '',
	'SMARTS': 'c1([SH1])aac([SH1])aa1',
	'redux': True,
},
{	
	'name': '1,2 - Aminophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c([OH1])c[NX3]([#1,#6])[#1,#6]',
	'redux': True,
},
{	
	'name': '1,3 - Aminophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c([OH1])[aR1]c[NX3]([#1,#6])[#1,#6]',
	'redux': True,
},
{	
	'name': '1,4 - Aminophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c1([OH1])aac([NX3]([#1,#6])[#1,#6])aa1',
	'redux': True,
},
{	
	'name': '1,2 - Aminothiophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c([SH1])c[NX3]([#1,#6])[#1,#6]',
	'redux': True,
},
{	
	'name': '1,3 - Aminothiophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c([SH1])[aR1]c[NX3]([#1,#6])[#1,#6]',
	'redux': True,
},
{	
	'name': '1,4 - Aminothiophenols',
	'description': 'R1, R2 = H, any carbon atom',
	'SMARTS': 'c1([SH1])aac([NX3]([#1,#6])[#1,#6])aa1',
	'redux': True,
},
{	
	'name': 'Alkynyl ethers',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX2]#[CX2][OX2][#6&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([OX2])=[OX1,SX1,NX2,C])]',
	'redux': True,
},
{	
	'name': 'Alkynyl thioethers',
	'description': 'R1, R2 = H, alkyl, aryl',
	'SMARTS': '[CX2]#[CX2][SX2][#6&!$([CX4]([OX2])([#7,O,S,F,Cl,Br,I,P]))&!$([CX3]([OX2])=[OX1,SX1,NX2,C])]',
	'redux': True,
},
{	
	'name': 'Carbonic acid esters',
	'description': 'R1 = H, alkyl, aryl; R2 = alkyl, aryl',
	'SMARTS': '[#1,#6&!$([CX3]=[OX1,SX1,NX2])][OX2][CX3](=[OX1])[OX2][#6&!$([CX3]=[OX1,SX1,NX2])]',
	'redux': True,
},
{	
	'name': 'Sulfuric acid esters',
	'description': 'R1 = alkyl, aryl; R2 = H, alkyl, aryl',
	'SMARTS': '[#6][OX2][Sv6X4](=[OX1])(=[OX1])[OH1,OX1-,$([OX2][#6])]',
	'redux': True,
},
{	
	'name': 'Phosphonic acid amides',
	'description': 'R1 = alkyl, aryl; R2, R3 = H, alkyl, aryl; Y = any O, N, Hal residue',
	'SMARTS': '[OX1]=[Pv5X4]([NX3])([#6])[!#6!$([SX2])]',
	'redux': True,
},
{	
	'name': 'Phosphonic acid halides',
	'description': 'R1 = alkyl, aryl; X = F, Cl, Br, I; Y = any O, N, Hal residue',
	'SMARTS': '[O,X1]=[Pv5X4]([F,Cl,Br,I])([#6])[!#6!$([SX2])]',
	'redux': True,
}
]
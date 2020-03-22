#!/usr/bin/env python3

import re
import itertools
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole

import rpy2
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

print('ssss')

mat_df  = pd.read_hdf('test.h5', key='df')
mat_df  = pd.read_csv('test.csv', delimiter=',', index_col=0)
mat_rdf = pandas2ri.py2ri(mat_df)
print(mat_rdf)


###############################################
from collections import defaultdict
from rdkit.Chem.rdchem import EditableMol
# Thanks to steeveslab-blog for example of how to edit RDKit molecules
# http://asteeves.github.io/blog/2015/01/14/editing-in-rdkit/
# Thanks to Andrew Dalke for the function name

from rdkit.Chem import rdDistGeom
m = Chem.MolFromSmiles('c1(N(C)C)ccc(N)cc1')
m2 = Chem.AddHs(m)
AllChem.EmbedMolecule(m2)
m_all = AllChem.EmbedMultipleConfs(m2, numConfs=10, numThreads=0)
len(m_all)

list(m_all)
rmslist = []
AllChem.AlignMolConformers(m2, RMSlist=rmslist)
len(rmslist)

rms = AllChem.GetConformerRMS(m2, 1, 9, prealigned=True)
rms
m2


def weld_r_groups(input_mol):
    # First pass loop over atoms and find the atoms with an AtomMapNum
    join_dict = defaultdict(list)
    for atom in input_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            join_dict[map_num].append(atom)

    # Second pass, transfer the atom maps to the neighbor atoms
    for idx, atom_list in join_dict.items():
        if len(atom_list) == 2:
            atm_1, atm_2 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()][0]
            nbr_1.SetAtomMapNum(idx)
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()][0]
            nbr_2.SetAtomMapNum(idx)

    # Nuke all of the dummy atoms
    new_mol = Chem.DeleteSubstructs(input_mol, Chem.MolFromSmarts('[#0]'))

    # Third pass - arrange the atoms with AtomMapNum, these will be connected
    bond_join_dict = defaultdict(list)
    for atom in new_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict[map_num].append(atom.GetIdx())

    # Make an editable molecule and add bonds between atoms with correspoing AtomMapNum
    em = EditableMol(new_mol)
    for idx, atom_list in bond_join_dict.items():
        if len(atom_list) == 2:
            start_atm, end_atm = atom_list
            em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.SINGLE)

    final_mol = em.GetMol()

    # remove the AtomMapNum values
    for atom in final_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    final_mol = Chem.RemoveHs(final_mol)

    return final_mol

###############################################
def anum(m, pat):
  return len(m.GetSubstructMatches(pat, maxMatches=m.GetNumAtoms()))

def Unsaturate(m):
  m     = Chem.AddHs(m)
  c_pat = Chem.MolFromSmarts('[#6]')             # carbon
  n_pat = Chem.MolFromSmarts('[#7]')             # nitrogen
  h_pat = Chem.MolFromSmarts('[#1,#9,#17,#35]')  #hydrogen and halogens
#  o_pat = Chem.MolFromSmarts('[#8]')             ## oxygen, no effect

  u = (2*anum(m,c_pat)+2 + anum(m,n_pat) - anum(m,h_pat))/2
  return int(u)

###############################################
test = list(range(0,5))
list(itertools.product(test, repeat=1))

ref = Chem.MolFromSmiles('C1NC=NC1')
ref
AllChem.Compute2DCoords(ref)  # storing 2D info into ref object
Chem.rdMolDescriptors.CalcMolFormula(ref)
Unsaturate(ref)
if re.search(r'\[\*:[0-9]\]', '[*:1]C([*:2])NC=NC([*:3])'): print('x')
num = re.findall(r'\[\*:(.*?)\]','[*:1]C([*:2])NC=NC([*:3])')

xxx = '[*:1]C([*:2])NC=NC([*:3])'
re.sub(r'\[\*:(.*?)\]', '%9'+str(num), '[*:1]C([*:2])NC=NC([*:3])')

r_group = ['CCC([*:1])', 'C(C)(C)([*:1])']
list_r  = [Chem.MolFromSmiles(m) for m in r_group]


s1 = Chem.MolFromSmiles('C1([*:1])NC=NC1')
s2 = Chem.MolFromSmiles('CCC([*:1])')
s3 = Chem.MolFromSmiles('C1(NCC)NC=NC1')
s3
Chem.CanonSmiles(s1+'.'+s2)

## Traditional way to add R-group to scaffod at specified locations
## %9{num} as flag for attachment site, instead of [*:{num}] as dummy
## %9{num} cannot be at beginning of SMARTS string, i.e. %91CCCCN=C failed
Chem.MolFromSmiles(Chem.CanonSmiles('C1%91NC=NC1'+'.'+'CCC%91'))
Chem.MolFromSmiles(Chem.CanonSmiles('C1%91%92NC=NC1%93'+'.'+'CCC%91.C%92.C(C)(C)%93'))
Chem.MolFromSmiles(Chem.CanonSmiles('C1%91CC2N=CNC12.CCCCCC%91'))
Chem.MolFromSmiles(Chem.CanonSmiles('C1%91CC2NC=NC2C1.%91C'))   # failed, %9 at start of string




##! This welding method by Patrick has issue with Geminal Carbon !##
## making it unusable for those cases, fallback to CanonSmiles %9 method
sr = Chem.MolFromSmiles('[*:1]C([*:2])NC=NC([*:3]).C(C)(C)C[*:1].C[*:2].C[*:3]')
AllChem.Compute2DCoords(sr)
sr
Chem.MolFromSmiles(Chem.MolToSmiles(weld_r_groups(sr)))

#############################################################
## single point
##  C1xxNC=NC1xx    # starting
scaff_s1 = ['C1([*:1])NC=NC1']
list_s1  = [Chem.MolFromSmiles(m) for m in scaff_s1]
## To align all structures to a reference structure, check if mol has matching substructure
## then compare to the reference's substructure 2D coords, update the coords of matched 
## substructure to the reference orientation
[AllChem.GenerateDepictionMatching2DStructure(m, ref) for m in list_s1]
Draw.MolsToGridImage(list_s1, molsPerRow=3,subImgSize=(125,125))

## double point
scaff_db = ['C1([*:1])NC=NC1([*:2])', 'C1([*:1])([*:2])NC=NC1']
list_db  = [Chem.MolFromSmiles(m) for m in scaff_db]
[AllChem.GenerateDepictionMatching2DStructure(m, ref) for m in list_db]
Draw.MolsToGridImage(list_db, molsPerRow=3)

## cyclopentyl and cyclohexyl
Chem.MolFromSmiles('C1CC2NC=NC2C1')
##  C1xxCxxC2xN=CNC12x          # cyclobutyl
scaff_cb = ['C1([*:1])CC2N=CNC12', 
            'C1([*:1])C([*:2])C2N=CNC12','C1([*:1])([*:2])CC2N=CNC12',
            'C1([*:1])([*:2])C([*:3])C2N=CNC12',
            'C1([*:1])([*:2)]C([*:3])([*:4])C2N=CNC12']
list_cb  = [Chem.MolFromSmiles(m) for m in scaff_cb]
[AllChem.GenerateDepictionMatching2DStructure(m, ref) for m in list_cb]
Draw.MolsToGridImage(list_cb, molsPerRow=3)

##  C1xxCxxC2xNC=NC2xC1xx       # cyclopentyl
scaff_cp = ['C1([*:1])CC2NC=NC2C1','C1C([*:1])C2NC=NC2C1',
            'C1([*:1])C([*:2])C2NC=NC2C1','C1C([*:1])C2NC=NC2C1[*:2]','C1([*:1])([*:2])CC2NC=NC2C1','C1C([*:1])([*:2])C2NC=NC2C1',
            'C1([*:1])C([*:2])C2NC=NC2C1([*:3])','C1([*:1])([*:2])C([*:3])C2NC=NC2C1','C1([*:1])C([*:2])([*:3])C2NC=NC2C1','C1C([*:1])([*:2])C2NC=NC2C1([*:3])',
            'C1([*:1])C([*:2])([*:3])C2NC=NC2C1([*:4])','C1([*:1])([*:2])C([*:3])C2NC=NC2C1([*:4])',]
list_cp  = [Chem.MolFromSmiles(m) for m in scaff_cp]
[AllChem.GenerateDepictionMatching2DStructure(m, ref) for m in list_cp]
Draw.MolsToGridImage(list_cp, molsPerRow=3)


##  C1xxCxxCxxC2xN=CNC2xC1xx    # cyclohexyl
scaff_ch = ['C1([*:1])CCC2N=CNC2C1','C1CCC2N=CNC2C1([*:1])',

            'C1([*:1])C([*:2])CC2N=CNC2C1','C1([*:1])CC([*:2])C2N=CNC2C1','C1([*:1])CCC2N=CNC2C1([*:2])','C1CC([*:1])C2N=CNC2C1([*:2])',
            'C1([*:1])([*:2])CCC2N=CNC2C1','C1CC([*:1])([*:2])C2N=CNC2C1',

            'C1([*:1])C([*:2])C([*:3])C2N=CNC2C1','C1([*:1])CC([*:2])C2N=CNC2C1([*:3])',
            'C1([*:1])([*:2])C([*:3])CC2N=CNC2C1','C1([*:1])([*:2])CC([*:3])C2N=CNC2C1',
            'C1([*:1])([*:2])CCC2N=CNC2C1([*:3])','C1([*:1])CC([*:2])([*:3])C2N=CNC2C1',
            'C1C([*:1])C([*:2])([*:3])C2N=CNC2C1','C1CC([*:1])([*:2])C2N=CNC2C1([*:3])',
            'C1([*:1])C([*:2])C([*:3])C2N=CNC2C1([*:4])'
]

scaff_ch = [
            'C1CCC2([*:1])N=CNC2C1',

            'C1([*:1])CCC2N=CNC2([*:2])C1','C1C([*:1])CC2N=CNC2([*:2])C1',
            'C1CC([*:1])C2N=CNC2([*:2])C1','C1CCC2N=CNC2([*:1])C1([*:2])','C1CCC2([*:1])N=CNC2([*:2])C1',

            'C1([*:1])CCC2N=CNC2([*:2])C1([*:3])','C1C([*:1])CC2N=CNC2([*:2])C1([*:3])',
            'C1CC([*:1])C2N=CNC2([*:2])C1([*:3])',
            'C1([*:1])C([*:3])CC2N=CNC2([*:2])C1','C1([*:1])CC([*:3])C2N=CNC2([*:2])C1',
            'C1([*:1])CCC2([*:3])N=CNC2([*:2])C1','C1CCC2([*:1])N=CNC2([*:2])C1([*:3])',

            'C1CCC2N=CNC2([*:2])C1([*:1])([*:3])','C1C([*:1])([*:3])CC2N=CNC2([*:2])C1',
            'C1CC([*:1])([*:3])C2N=CNC2([*:2])C1','C1([*:1])([*:3])CCC2N=CNC2([*:2])C1',
]
scaff_ch = [
            'C1([*:1])C([*:2])CC2([*:3])N=CNC2([*:4])C1',
            'C1([*:1])CC([*:2])C2([*:3])N=CNC2([*:4])C1','C1([*:1])([*:2])CCC2([*:3])N=CNC2([*:4])C1',
            'C1CC([*:1])([*:2])C2([*:3])N=CNC2([*:4])C1',

            'C1([*:1])([*:3])C([*:2])([*:4])CC2N=CNC2C1',
            'C1([*:1])([*:2])CC([*:3])([*:4])C2N=CNC2C1','C1([*:1])([*:2])CCC2N=CNC2C1([*:3])([*:4])',
            'C1CC([*:1])([*:2])C2N=CNC2C1([*:3])([*:4])',

            'C1([*:1])CC([*:3])C2N=CNC2([*:2])C1([*:4])','C1CC([*:3])C2N=CNC2([*:2])C1([*:1])([*:4])',
            'C1C([*:1])C([*:3])C2N=CNC2([*:2])C1([*:4])','C1CC([*:1])([*:3])C2N=CNC2([*:2])C1([*:4])',
            ]

list_ch  = [Chem.MolFromSmiles(m) for m in scaff_ch]
[AllChem.GenerateDepictionMatching2DStructure(m, ref) for m in list_ch]
Draw.MolsToGridImage(list_ch, molsPerRow=4, subImgSize=(125,125))

ch_df = pd.DataFrame(scaff_ch, columns=['core'])
ch_df.attach[0]
ch_df['attach'] = ch_df.core.map(lambda m: m.count('*'))


## generate all permutations of attachment points on scaffold with 
## attachment designation
all = []
s = list('C1xxCxxCxxC2xN=CNC2xC1xx')
yyy = [list(map(int,x)) for x in list(itertools.product('01', repeat=10))]
xxx = [y for y in yyy if sum(y) <= 4]
len(xxx)
for xx in xxx:
  mol = s
  for i, x in enumerate(xx):
    for j, m in enumerate(mol):
      if m =='x':
        if x:
          print(i,j)
          mol = ''.join(mol[:j])+'([*])'+''.join(mol[j+1:])
          break
        else:
          mol = ''.join(mol[:j])+''.join(mol[j+1:])
          break
  all.append(mol)
  print(mol)
zzz = set([Chem.MolToSmiles(Chem.MolFromSmiles(m)) for m in set(all)])
len(zzz)
len(set(zzz))
Draw.MolsToGridImage(zzz, molsPerRow=4, subImgSize=(125,125))

for m in set(all):
  Chem.MolFromSmiles(m)

len(all)
for a in all:
  print(a)

a = 'CC12CC'
b = 'CC[C@H]CC'
c = '[*][C@@H](N)(O)CC'
d = 'C[*]12CCCC12'
e = 'CCC[CH2+]CCC'
if re.search('[1-9]+',a): print('True')
re.sub('\[\*\]([1-9]+)', 'a[*]',d)
if re.search('\[CH2\]|\[CH\]|\[CH2\+\]', e): print('True')

z= [('x',0,'*C\=/CCC1C2'),('x',0,'[*]NCCC'),('x',0,'N*CCC'),('x',0,'[*][C@@H]1(N)(O)C1C')  ]

collect = []
for m in z:
  if not re.search(r'\[\*\]', m[2]):
    n = re.sub('\*', '[*]', m[2])
    print(n)
    m = ( m[0], m[1], n )   # save as tuple for set() later
    
  # if [*] is at beginning, [*]CCCC
  if re.search(r'^\[\*\]', m[2]):
    head = re.findall(r'^\[\*\](.?)',  m[2])[0]
    # if what follows [*] is [???], ie [C@H]
    if re.search(r'\[', head):
      head = re.findall(r'\[\*\](\[.+?\])', m[2])[0]
      temp = re.sub('\[\*\]\[.+\]', head+'[*]', m[2])
      m = (m[0], m[1],temp )
    else:
      temp = re.sub('\[\*\]([A-Za-z0-9])', head+'[*]', m[2])
      m = (m[0], m[1], temp)
    
  if re.search('\[\*\]([1-9]+)', m[2]):
    head = re.findall(r'\[\*\]([1-9]+)', m[2])[0]
    temp = re.sub('\[\*\]([1-9]+)', head+'[*]', m[2])
    m = (m[0], m[1], temp)

  collect.append(m)
  print(collect)





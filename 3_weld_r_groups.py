#!/usr/bin/env python3

import sys,re
import bz2
import pickle
import itertools
import numpy as np
import pandas as pd

from rdkit import rdBase
rdBase.rdkitVersion

from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions

from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser
from collections import defaultdict
chunk = 500   # chunk size for mutliprocessing

##########################################################################

def main():
  args  = UserInput()

  if args.unsat_max:
    unsat_max = int(args.unsat_max)
  else:
    unsat_max = None
  if args.unsat_min:
    unsat_min = int(args.unsat_min)
  else:
    unsat_min = None
  if args.max_c:
    max_c = int(args.max_c)
  else:
    max_c   = 11            # hard limit of carbon atom for Mol-X
  if args.allow_c:
    allow_c = int(args.allow_c)
  else:
    allow_c = 8             # number of extra carbon allowed for R-group


  ## the attachment points are annotated with '%9' instead of [*], i.e.
  ## '%91' and '%93' is equivalent to describing attachment points [*:1]
  ## and [*:3] in SMARTS format. However since Chem.CanonSmiles() function
  ## only takes '%9' flag, all input [*:{sum}] will be replaced with '%9'
  ##
  ## Unlike SMARTS annotation where () bracket is needed fro branching pt,
  ## e.g. C(C)([*])C([*][*])C, '%9' method does not need (), i.e. C(C)%91C%92%93C
  mpi = multiprocessing.Pool()

################  Read data into dataframe  ##################

  # Read in plain text of cores, make sure attachment points are in 'x'
  # format, then account for number of attachment point
  # remove duplicates, cannot remove tautomers at this stage
  with open(args.tmpl_file, 'r') as fi:
    Templates = [l.strip() for l in fi]
  Temp_1 = [GenCoreFromString(inp, max_c) for inp in Templates]
  Smiles = [item for sublist in Temp_1 for item in sublist] # flatten
  print(' # Initial core generated: '+str(len(Smiles)))

  core_df = pd.DataFrame(Smiles, columns=['smiles','canon']).drop_duplicates('smiles').drop_duplicates('canon')
  core_df['core']   = core_df.smiles.map(lambda m: FlagRegex(m))
  core_df['attach'] = core_df.core.map(lambda m: m.count('%9'))
  core_df['maxatm'] = core_df.core.map(lambda m: max_c - m.count('C'))
  print(core_df[:10])
  print(' # Nonredundant core: '+str(len(core_df)))


########  For each core, exhaustive permutations of frag decoration  #########

  if args.raw is None:
    # put R-group into dataframe [r_group, smiles, atom_num]
    r_df = pd.read_csv(args.r_file, sep='\s+', header=0, comment='#')
    print(r_df)
    print(' # Nonredundant frag: '+str(len(r_df)))

    Mols = []
    ## iterate every row if core dataframe, add R-group(s) to attachment 
    ## point(s) on core fragments
    for idx, row in core_df.iterrows():
      print('# {0} of {1}: {2}'.format(idx+1, len(core_df), row.core ))

      # Single attachment, add R-group directly to core
      if row.attach == 1:
        # limit R-group to those won't exceed row.maxatom
        rg_mols = r_df[ r_df.atom_num == row.maxatm ].r_group
        for rnum, rg in enumerate(rg_mols):
          name = 'mol_{0}_{1}_{2}'.format(idx+1, row.attach, rnum+1)
          Mols.append([ row.core+'.'+re.sub(r'\[\*\]', '%91', rg ), name ])

      # Multiple attachment, generate R-group combos before add to core
      else:
        # limit R-group number to those won't exceed upper boundary
        rg_mols = r_df[r_df.atom_num <= row.maxatm-row.attach+1].r_group

      # generate list of r-group combo with specified Carbon number
      rg_list = RGroupProduct(rg_mols, row.attach, row.maxatm)
      for rnum, rg in enumerate(rg_list):
        name = 'mol_{0}_{1}_{2}'.format(idx+1, row.attach, rnum+1)
        Mols.append([ row.core + rg, name ])
    l = len(Mols)

    CombineFrag = CombineFragments()
    ## Combine core/R-group fragments to generate molecule
#    Tmp = [x for x in tqdm(mpi.map(CombineFrag, Mols, chunk), total=l)]
    Tmp = [CombineFrag(m) for m in Mols]
    Analogs = [m for m in Tmp if m is not None]


    ## Save/load Analog intermediate reuslts to pickle
    with bz2.BZ2File(args.outpref+'.pickle.bz2', 'w') as fpo:
      pickle.dump(Analogs, fpo, protocol=pickle.HIGHEST_PROTOCOL)

  else:
    with bz2.open(args.raw, 'rb') as fpi:
      Analogs = pickle.load(fpi)


###############  proprocess generated initial combinations #################

  ## Remove duplicate mols and limit deg of unsaturation
  print('## Removing duplicated SMILES from {0}...'.format(len(Analogs))) 
  Cols = ['smiles', 'name', 'unsat', 'formula', 'molwt']
  a_df = pd.DataFrame(Analogs, columns=Cols).drop_duplicates('smiles')

  if unsat_max is not None or unsat_min is not None:
    print('## Limiting degree of unsaturation from {0}...'.format(len(a_df)))
    if unsat_max is not None:
      a_df = a_df[ a_df.unsat < unsat_max ]
    if unsat_min is not None:
      a_df = a_df[ a_df.unsat > unsat_min ]

  ## Generate unique stereoisomers (double bonds, chiral centers, no meso-)
  print('## Generating unique stereoisomers (no meso-isomers) from {0}...'.format(len(a_df)))
  Mols     = a_df.values.tolist()
  l        = len(Mols)
  Temp_1   = [x for x in tqdm(mpi.imap(Stereoisomer, Mols, 5),total=l)]
#  Temp_1 = [Stereoisomer(mol) for mol in Mols]
  Sto_Mols = [item for sublist in Temp_1 for item in sublist]   # flatten
  print('## Summarize results...')
  t_df     = pd.DataFrame(Sto_Mols, columns=Cols).drop_duplicates('smiles')
  final_df = t_df.sort_values(['unsat','smiles'])
  mpi.close()
  mpi.join()


#################  Summarize and output results  ################

  print('\n  ## Core Generated:    {0:10d} -- Unique: {1}\n'.format(
           len(Smiles), len(core_df)))
  print('\n  ## Total Permutation: {0:10d} -- Unique: {1}\n'.format(
           len(Analogs), len(a_df) ))
  print('\n  ## Total Steroisomer: {0:10d} -- Unique: {1}\n'.format(
           len(Sto_Mols), len(final_df) ))

  final_df.to_csv(args.outpref+'.txt', sep='\t', header=True, index=False)
  final_df.to_csv(args.outpref+'.smi', sep='\t', header=True, index=False,
                  columns=['smiles', 'name'])


##########################################################################
##########################################################################
## Put the attachment Site designation onto the R-groups by replacing [*]
## SMARTS flag with %9{num} flag for Chem.CanonSmiles()
def RGroupProduct( mols, attach, maxatm ):

  class JoinRFrag(object):
    def __init__( self, maxatm = 0 ):
      self.maxatm = maxatm    # maximum allowed (needed) C atom

    def __call__( self, combo ):
      return self.run_combine( combo )

    # substitute the dummy from R-group with %9 flag then combine strings
    # only return strings matching maximum allowed C atom
    def run_combine( self, combo ):
      linked = ''
      for idx, m in enumerate(combo):
        linked += '.'+re.sub(r'\[\*\]', '%9'+str(idx+1), m)
      if linked.count('C') == self.maxatm:
        return linked
      else:
        return None

  ## Create an exhaustive list of R-group permutations, including self.self
  ## then in JoinRFrag reduce the list to those with maximum allowed number 
  ## of carbon atom 'allow_c'
  Mol = list(itertools.product(mols, repeat=attach))
  mpi = multiprocessing.Pool()
  frg = JoinRFrag(maxatm=maxatm)
#  Tmp = [x for x in tqdm(mpi.imap(frg, Mol, chunk),total=len(Mol))]
  Tmp = [frg(m) for m in Mol]
  mpi.close()
  mpi.join()

  return [r_group for r_group in Tmp if r_group is not None]
  

##########################################################################
## Combine Core/R-group fragments into one single molecule and get molecule
## properties. Welding fragments together using Chem.CanonSmiles() works.
## Welded molecule is standardized, tautomerized, and canonicalized so that
## all can be compared directly for duplication
class CombineFragments(object):
  def __init__( self ):
    self.c_pat = Chem.MolFromSmarts('[#6,#14]')         # carbon, silicon
    self.n_pat = Chem.MolFromSmarts('[#7,#15]')         # nitrogen, Phosphorus
    self.h_pat = Chem.MolFromSmarts('[#1,#9,#17,#35]')  # hydrogen and halogens

  def __call__( self, inp ):
    return self.combine_frags(inp)

  def combine_frags( self, inp ):
    fragments, name = inp   # fragments = 'C1%91C%92NC=N1.C%91C.C1%92CC1'

  #  tau_mol = _weld_mol(fragments)   # failed SMARTS method
    tau_mol = _comb_mol(fragments)
    tau_smi = Chem.MolToSmiles(tau_mol)
    dgunsat = DegOfUnsaturation(tau_mol, self.c_pat, self.n_pat, self.h_pat)
    formula = Chem.rdMolDescriptors.CalcMolFormula(tau_mol)
    molwt   = Chem.rdMolDescriptors.CalcExactMolWt(tau_mol)

    return [ tau_smi, name, dgunsat, formula, molwt ]


##############################
# welding with CanonSmiles() %9(num) flag, a stable and safe method
# it takes a core and multiple fragments, with attachment points flagged by
# %9(num), in sequential format, i.e. frags = 'C1%91C%92NC=N1.C%91C.C1%92CC1'
# Chem.CanonSmiles joins them together
# Standardize, tautomerize, canonicalize the structure for later comparison
def _comb_mol( frags ):
  print(frags)
  std_smi = MolStandardize.standardize_smiles(Chem.CanonSmiles(frags))
  tau_smi = MolStandardize.canonicalize_tautomer_smiles(std_smi)
  tau_mol = Chem.MolFromSmiles(tau_smi)
  return tau_mol

# welding with SMARTS input [*:x], but failed geminal input
def _weld_mol( frags ):
  combine = Chem.MolToSmiles( Weld_R_Groups( Chem.MolFromSmiles(frags) ))
  std_smi = MolStandardize.standardize_smiles(combine)
  tau_smi = MolStandardize.canonicalize_tautomer_smiles(std_smi)
  tau_mol = Chem.MolFromSmiles(tau_smi)
  return tau_mol


##########################################################################
## Generate stereoisomers then tautomerize molecules to make sure
## equivalent stereoisomers have the same Smiles output
def Stereoisomer( inp ):
  smiles, name, dgunsat, formula, molwt = inp

  sto_opt = StereoEnumerationOptions(tryEmbedding=False, unique=True)

  mol     = Chem.MolFromSmiles(smiles)
  ism_mol = tuple( EnumerateStereoisomers(mol, options=sto_opt) )
  ism_smi = [Chem.MolToSmiles(s, isomericSmiles=True) for s in ism_mol]
  tau_smi = [MolStandardize.canonicalize_tautomer_smiles(m) for m in ism_smi]
  
  Out = []
  for idx, smi in enumerate(tau_smi):
    Out.append([smi, name+'_'+str(idx+1), dgunsat, formula, molwt])

  return Out


###############################################
##  get number of specified atom (SMARTS pattern)
def anum( m, pat ):
  return len(m.GetSubstructMatches(pat, maxMatches=m.GetNumAtoms()))

## calculate degree of unsaturation of input molecule
## oxygen/phosphous do not affect DoU calculation
def DegOfUnsaturation( m ):
  c_pat = Chem.MolFromSmarts('[#6]')             # carbon
  n_pat = Chem.MolFromSmarts('[#7]')             # nitrogen
  h_pat = Chem.MolFromSmarts('[#1,#9,#17,#35]')  # hydrogen and halogens

  m = Chem.AddHs(m)
  u = (2*anum(m,c_pat)+2 + anum(m,n_pat) - anum(m,h_pat))/2
  return int(u)


##########################################################################
## using 'x' as flag for all possible attachment points on any atom as initial
## read-in, generate SMARTS format with maximum number of attach points the
## scaffold allows with the specific limit of R-group/core carbon numbers
## example of 'x' flagged scaffold: mol = 'C1xxCxxCxxC2xN=CNC2xC1xx'

def GenCoreFromString( mol, max_c ):

  repeat = mol.count('x')
  c_num  = mol.count('C')
  Temp_1 = [x for x in list(itertools.product('01', repeat=repeat))]
  Temp_2 = [list(map(int, t)) for t in Temp_1]    # convert str to int

  # limit attachment points to maximum the scaffold allows
  Product= [t for t in Temp_2 if sum(t) <= (max_c - c_num)]

  Cores = []
  for combo in Product:   # list of [1,0,0,0,1,1,1,0]
    m = mol
    n = mol
    for sub in combo:
      # molecule for final use, [*] has no ()
      for idx, atom in enumerate(n):
        if atom == 'x':
          if sub:
            n = ''.join(n[:idx])+'[*]'+''.join(n[idx+1:])
            break
          else:
            n = ''.join(n[:idx])+''.join(n[idx+1:])
            break
      # For canonicalization use, [*] has (), not final use
      for idx, atom in enumerate(m):
        if atom == 'x':
          if sub:
            m = ''.join(m[:idx])+'([*])'+''.join(m[idx+1:]) # for canon use
            break
          else:
            m = ''.join(m[:idx])+''.join(m[idx+1:])
            break

    cleaned = CleanFragBranch(n)    # kinda unnecessary with the current setup
    if re.search(r'\*', cleaned):   # save molecules with attach point(s)
      ## get a canonicalized version of the molecule for comparison
      std_smi = MolStandardize.standardize_smiles(m)
      tau_smi = MolStandardize.canonicalize_tautomer_smiles(std_smi)
      Cores.append([cleaned, tau_smi])

  return Cores


##########################################################################
## Check if input core fragment has attachment point and in correct format
## if found ([*]), ignore ()
## Convert dummy attachment points 'x' into numbers
def FlagRegex(inp):
  ## if attachment points are in SMARTS format and numbered [*:1]
  if re.search(r'\*:[0-9]', inp):
    inp = re.sub(r'\(\[\*:', '%9', inp)   # C[*:1]([*:2])NC=N -> C[*:1]%92NC=N
    inp = re.sub(r'\[\*:', '%9', inp)     # C[*:1]([*:2])NC=N -> C%91([*:2])NC=N
    inp = re.sub(r'\]\)', '', inp)
    inp = re.sub(r'\]', '', inp)
    return inp
  ## if attachment points are only in SMARTS format
  elif re.search(r'\[\*\]', inp):
    num = list(range(1,inp.count('*')+1))
    inp = re.sub(r'\(\[\*\]\)', '%9x', inp) # C[*]([*])NC=N -> C[*]%9xNC=N
    inp = re.sub(r'\[\*\]', '%9x', inp)     # C[*]([*])NC=N -> C%9x([*])NC=N
    for i in range(len(inp)):
      if inp[i] == 'x':
        inp = inp[:i]+ str(num.pop(0)) +inp[i+1:]
    return inp
  ## if attachment points are already in %9 format and numbered
  elif re.search(r'\%9', inp):
    return inp
  else:
    sys.exit('  ## ERROR: No branch point (%9) found in core fragment: '+inp)


##########################################################################
## clean up cases of no bracket: *C=C -> [*]C=C
## canonical smiles places [*] at beginning of string but this will fail
## Chem.CanonSmiles()'s %9 flag, need to place something (carbon) before
## [*] flag here. Also need to deal with branch and special atoms, e.g.
## C12, [C@@H], [C@H]12 
## for [*] in smiles format, they need to have () to designate branch
## but for %9, no () is needed

def CleanFragBranch( frag ):
  ## if found * without bracket, add [] to become [*]
  if not re.search(r'\[\*\]', frag):
    frag = re.sub('\*', '[*]', frag)


  ## if [*] at beginning, [*]CCCC, move it to *after* first atom
  if re.search(r'^\[\*\]', frag):
    head = re.findall(r'^\[\*\](.?)', frag)[0]

    # if follows [*] is [???], i.e. [C@H], add () and move again
    if re.search(r'\[', head):
      head = re.findall(r'^\[\*\](\[.+?\])', frag)[0]
      print(head)
      print(frag)
      frag = re.sub(r'^\[\*\](\[.+?\])', head+'[*]', frag)
      print(frag)
    else:
      frag = re.sub(r'\[\*\]([A-Za-z0-9])', head+'[*]', frag)

    # if follows [*] is double '/X' '/[X]', move it to *after* first atom
    if re.search(r'\/', head):
      head = re.findall(r'^\[\*\](\/\[.+?\])', frag)[0]
      print(head)
      print(frag)
      frag = re.sub(r'^\[\*\](\/\[.+?\])', head+'[*]', frag)
      print(frag)
    else:
      frag = re.sub(r'\[\*\](\/[A-Za-z0-9])', head+'[*]', frag)


  ## if follows [*] is number, i.e. C[*]1+, move again
  if re.search('\[\*\]([1-9]+)', frag):
    head = re.findall(r'\[\*\]([1-9]+)', frag)[0]
    frag = re.sub('\[\*\]([1-9]+)', head+'[*]', frag)


  ## if follows [*] is double bond, i.e. '/''\', move again
  ## there may be multiple double bonds, need looping
  if re.search(r'\[\*\]\\', frag):
    heads = re.findall(r'\[\*\](\\.?)', frag)
    for head in heads:
      frag = re.sub(r'\[\*\](\\.?)', head+'[*]', frag)
  if re.search('\[\*\]\/', frag):
    heads = re.findall(r'\[\*\](\/.?)', frag)
    for head in heads:
      frag = re.sub('\[\*\](\/.?)', head+'[*]', frag)


  return frag


##########################################################################
###############################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-templ', dest='tmpl_file', required=True,
                 help='Core Scaffold SMILES with attachment points marked by "x"\n (eg: CxxCx=C)')
  p.add_argument('-out', dest='outpref', required=True,
                 help='Output prefix')

  p.add_argument('-r', dest='r_file', required=False,
                 help='CSV file of R groups with attachment point marked by "[ * ]" generated by 2_parse_mmpdb_frag.py')
  p.add_argument('-raw', dest='raw', required=False,
                 help='Optional: Read-in pre-generated analog intermediate file, i.e. pickle file of R groups (def: None)')

  p.add_argument('-unsat_min', dest='unsat_min', required=False,
                 help='Optional: Remove molecule with deg unsaturation less than this_(def: None)')
  p.add_argument('-unsat_max', dest='unsat_max', required=False,
                 help='Optional: Remove molecule with deg unsaturation larger than this (def: None)')

  p.add_argument('-max_c', dest='max_c', required=False,
                 help='Optional: Maximum carbon number for final product (def: 11)')
  p.add_argument('-allow_c', dest='allow_c', required=False,
                 help='Optional: Number of carbon allowed for R-group for welding def: 8)')


  args=p.parse_args()
  return args

##########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
#   Peter M.U. Ung @ MSSM/Yale
#
#   v1.0   19.05.12
#   v2.0   19.12.12  new cases of string naming
#
#   Create combinatorial analog library of a core scaffold using a fragment 
#   library While working on the generation, intermediates from the "Combine
#   Core/R-group to generate molecule" step is saved into a pickle.bz2 in case
#   failure in later steps.
#
#   1) create all permutations of cores with different number of branch pts
#       limited by predefined number of allowed R-group carbon atom
#   2) create all permutations of different number of all R-groups, limited
#       by predefined number of allowed R-group carbon atom
#   3) Combine core scaffolds (with designated branch points) and R-groups 
#       (with branch points) into one molecule
#   4) remove molecules with degree of unsaturation not fitting criteria
#   5) remove duplicated SMILES (strings are tautomerized and canonicalized)
#   6) generate all possible stereoisomers (E/Z, diastereomer) and remove 
#       meso-isomers (SMILES are tautomerized and canonicalized)
#   7) report results
#
#
##########################################################################
##! This welding method uses SMARTS [*] annotation directly !##
## works in most cases but FAILED for germinal attachment point, 
## i.e. CC([*:1])([*:2])CCC. Instead the classic Chem.CanonSmiles()
## method works for germinal attachment. Decided not to use this one

## https://sourceforge.net/p/rdkit/mailman/message/36294557/
# Thanks to steeveslab-blog for example of how to edit RDKit molecules
# http://asteeves.github.io/blog/2015/01/14/editing-in-rdkit/
# Thanks to Andrew Dalke for the function name
def Weld_R_Groups(input_mol):
  from rdkit import Geometry
  from rdkit.Chem import Draw
  from rdkit.Chem.rdchem import EditableMol
  from rdkit.Chem.Draw import IPythonConsole

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

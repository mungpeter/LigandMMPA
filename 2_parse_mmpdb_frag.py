#!/usr/bin/env python3

import sys
import re
import gc
import time
import json
import pickle
import gzip,bz2
import itertools
import pandas as pd
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser
chunk = 500     # mpi chunk size

from rdkit import Chem
from rdkit import rdBase

pd.__version__
rdBase.rdkitVersion

##########################################################################
def main():

  args = UserInput()
  if args.size:
    size  = int(args.size)
  else:
    size  = 999
  if args.regex:
    regex = args.regex
  else:
    regex = 'Se|Te|B'
  if args.unsat:
    unsat = int(args.unsat)
  else:
    unsat = 999

#######################
  with open(args.f_list, 'r') as fi:
    infile = list(filter(None, (l.rstrip() for l in fi)))
#  print(infile)

  Temp = [ ProcessFrag(fi, size, regex, unsat) for fi in infile ]
  r_df = pd.DataFrame(list( set(tuple(itertools.chain.from_iterable(Temp)))),
            columns=['r_group','smiles','atom_num'] ).sort_values(['atom_num'])

  print('## number of unique frag: ', len(r_df))
  r_df.to_pickle(args.outpref+'.pickle.bz2')
  r_df.to_csv(args.outpref+'.txt.bz2', sep='\t', header=True, index=False, compression='bz2')
  r_df.to_csv(args.outpref+'.smi.bz2', sep='\t', header=True, index=False, columns=['r_group','smiles'], compression='bz2')


######################################################################
######################################################################
## A class to read in and process the JSON-formatted fragment data, see end of
## this script for details
class read_in_json(object):

  def __init__( self, size='', regex='', unsat='' ):
    self.size  = size
    self.regex = regex
    self.unsat = unsat
    self.c_pat = Chem.MolFromSmarts('[#6,#14]')         # carbon, silicon
    self.n_pat = Chem.MolFromSmarts('[#7,#15]')         # nitrogen, phosphorus
    self.h_pat = Chem.MolFromSmarts('[#1,#9,#17,#35]')  # hydrogen and halogens

  def __call__( self, x ):
    return self._load(x)
  
  def _load( self, x ):
    y = json.loads(x)
    z = []

    ## Save into list of list: [ ['frag with [*]', 'smiles of frag', 'HA number' ], [], ...]
    if not re.search('IGNORE', y[0]):
      print('smiles: ', y[2])
      try:
        for i in y[5]:
          print(i)
          # if json[5][i][8] has separate frags, [9] will be 'null' and return None, 
          # also remove any radical and cations
          if i[9] is not None and not re.search(r'\[CH\]|\[CH2\]|\[CH2\+\]', i[8]):
            ## exclude fragment if match these regex pattern or too many unsaturations
            if i[6] <= self.size and not re.search(self.regex, i[9]):
              dou = DegOfUnsaturation(i[9], self.c_pat, self.n_pat, self.h_pat)
              if dou <= self.unsat:
                z.append( (i[8], i[9], i[6]) ) # save as tuple for set() later

        ## clean up bad cases of SMILES branch point:
        ## 1) no bracket: *C=C -> [*]C=C, 2) [*] at front of string [*]CCCCCCC,
        ## 3) [*] when in front of cyclic number C[*]12 or special atom [*][C@H]
        ## 4) no need for () with [*] as Chem.CanonSmiles()'s %9 flag will fail
        ## resave the cleaned frag with rest of info in 'tuple' for set() later
        Cleaned_frag = [ (CleanFragBranch(m[0]), m[1], m[2]) for m in z ]
        return Cleaned_frag

      except IndexError:
#        print(y)
        print('x')

########################################################
## process the json formatted fragment library and select fragments with criteria
def ProcessFrag( f, size, regex, unsat ):

  with file_handle(f) as fi:
#    print(f)
    Tmp = [l.rstrip() for l in fi]
  Frags = Tmp[10:]    # remove first 10 lines, which are just file formatting
  del Tmp

  mpi = multiprocessing.Pool(processes=multiprocessing.cpu_count())
  frag = read_in_json(size=size, regex=regex, unsat=unsat)

  start = time.perf_counter()
  Tmp = [frag(x) for x in Frags]
#  Tmp = [frag(x) for x in tqdm(Frags, total=len(Frags))]
#  Tmp = [x for x in tqdm(mpi.imap(frag, Frags, chunk),total=len(Frags))]
  mpi.close()
  mpi.join()
  end = time.perf_counter()
  print((end-start), 'ms')
  Tmp2 = [frag for frag in Tmp if frag is not None]
  Tmp3 = list(itertools.chain.from_iterable(Tmp2)) # flatten nested lists
  FgSele = list(set(tuple(Tmp3)))   # select unique smiles

  del Tmp
  del Tmp2
  del Tmp3
  del Frags
  gc.collect()  # collect "gabage" to free up memory

  return FgSele


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


  ## if [*] at beginning, [*]CCCC, move it after first atom
  if re.search(r'^\[\*\]', frag):
    head = re.findall(r'^\[\*\](.?)', frag)[0]

    # if follows [*] is [???], i.e. [C@H], always add () and move again
    if re.search(r'\[', head):
      head = re.findall(r'\[\*\](\[.+?\])', frag)[0]
      frag = re.sub(r'\[\*\]\[.+\]', head+'[*]', frag) 
    else:
      frag = re.sub(r'\[\*\]([A-Za-z0-9])', head+'[*]', frag)

    # if follows [*] is double '/X' '/[X]', move it to *after* first atom
    if re.search(r'\/', head):
      head = re.findall(r'^\[\*\](\/\[.+?\])', frag)[0]
      print(head)
#      print(frag)
      frag = re.sub(r'^\[\*\](\/\[.+?\])', head+'[*]', frag)
      print(frag)
    else:
      frag = re.sub(r'\[\*\](\/[A-Za-z0-9])', head+'[*]', frag)


  ## if follows [*] is number, i.e. 1+, move again
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

  print(frag)
  return frag


###############################################
##  get number of specified atom (SMARTS pattern)
def anum( m, pat ):
  return len(m.GetSubstructMatches(pat, maxMatches=m.GetNumAtoms()))

## calculate degree of unsaturation of input molecule
## oxygen/phosphous do not affect DoU calculation
def DegOfUnsaturation( mol, c_pat, n_pat, h_pat ):

  mol = Chem.AddHs(Chem.MolFromSmiles(mol))
  dou = (2*anum(mol,c_pat)+2 + anum(mol,n_pat) - anum(mol,h_pat))/2

  return int(dou)


##########################################################################
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    return gzip.open(file_name, 'r')
  elif re.search(r'.bz2$', file_name):
    return bz2.BZ2File(file_name, 'r')
  else:
    return open(file_name, 'r')

##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-list', dest='f_list', required=True,
                  help='List of fragmentation smiles files')

  p.add_argument('-out', dest='outpref', required=True,
                  help='Output prefix')

  p.add_argument('-size', dest='size', required=False,
                  help='Maximum fragment heavy atom count allowed (def: 999)')
  p.add_argument('-unsat', dest='unsat', required=False,
                  help='Maximum Degree of Unsaturation allowed (def: 999)')
  p.add_argument('-regex', dest='regex', required=False,
                  help='Regular Expression of atomtypes to be excluded\n(def: "c|n|s|N|S|O|P|F|Cl|Br|I|Se|Te|B")')
  args=p.parse_args()
  return args


##########################################################################
if __name__ == '__main__':
  main()


##########################################################################
#
#   Peter M.U. Ung @ MSSM/Yale
#
#   v1.0  19.05.10 
#   v2.0  19.12.12  remove some setting limitations
#
#   Parse the JSON results from mmpDB's fragmentation step originally for 
#   Matched Molecular Pair Analysis. 
#   Here the results are collected to generate a list of SMILES strings with
#   attachment point designated as [*].
#
#   [*] are placed at the beginning of the string by mmpDB, may create problems
#   for Chem.CanonSmiles() function when combining R-groups to core molecule.
#   Place the [*] flag behind one atom:
#     [*]CCCC          --> C([*])CCC
#     [*]C1CC1         --> C1([*])CC1
#     [*][C@H]12CC1C2  --> [C@H]12([*])CC1C2
#
##########################################################################
# the json output of mmpdb is a nest of lists with variable list size, making
# standard json readin not feasible. The read_in_json class here manually parse
# the nested lists format, examples:
#
# ["RECORD", "phenol", "Oc1ccccc1", 7, "Oc1ccccc1",
#   [ [1, "N", 1, "1", "[*]O", "0", 6, "1", "[*]c1ccccc1", "c1ccccc1"],
#     [1, "N", 6, "1", "[*]c1ccccc1", "0", 1, "1", "[*]O", "O"]             ]]
#
# ["RECORD", "catechol", "Oc1ccccc1O", 8, "Oc1ccccc1O",
#   [ [1, "N", 1, "1", "[*]O", "0", 7, "1", "[*]c1ccccc1O", "Oc1ccccc1"],
#     [1, "N", 7, "1", "[*]c1ccccc1O", "0", 1, "1", "[*]O", "O"],
#     [2, "N", 6, "11", "[*]c1ccccc1[*]", "01", 2, "11", "[*]O.[*]O", null] ]]
#
# ["RECORD", "2-aminophenol", "Oc1ccccc1N", 8, "Nc1ccccc1O",
#   [ [1, "N", 1, "1", "[*]N", "0", 7, "1", "[*]c1ccccc1O", "Oc1ccccc1"],
#     [1, "N", 1, "1", "[*]O", "0", 7, "1", "[*]c1ccccc1N", "Nc1ccccc1"],
#     [1, "N", 7, "1", "[*]c1ccccc1N", "0", 1, "1", "[*]O", "O"],
#     [1, "N", 7, "1", "[*]c1ccccc1O", "0", 1, "1", "[*]N", "N"],
#     [2, "N", 6, "11", "[*]c1ccccc1[*]", "01", 2, "12", "[*]N.[*]O", null] ]]
#
# json[0] - record of fragments from one molecule
# json[1] - molecule name
# json[2] - molecule SMILES
# json[3] - number of atoms in full molecule
# json[4] - ?? canonical smiles ??
# json[5] - list of fragments
# json[5][i] - individual fragment
# json[5][i][0] - number of cuts in fragment
# json[5][i][1] - enumerate label (irrelevant)
# json[5][i][2] - (variable) number of heavy atom in frag
# json[5][i][3] - (variable) symmetry class
# json[5][i][4] - (variable) smiles of frag with attachment point [*]
# json[5][i][5] - order of frag
# json[5][i][6] - (constant) number of heavy atom in frag
# json[5][i][7] - (constant) symmetry class
# json[5][i][8] - (constant) smiles of frag with attachment point [*]
# json[5][i][9] - (constant) smiles of frag with hydrogen
#
# In usage, the (constant) section of json is used, json[5][i][9] as test if frag
# has only 1 attachment point, and use json[5][i][8] as fragment
#

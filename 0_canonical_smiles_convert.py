#!/usr/bin/env python3

## convert obabel canonical smiles into rdkit canonical smiles

import sys,re,gc
import gzip,bz2
from rdkit import Chem

def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'r')
  elif re.search(r'.bz2$', file_name):
    print(file_name)
    handle = bz2.BZ2File(file_name, 'r')
  else:
    handle = file_name

  return handle

def rdkit_open(File_Tuple):

  List = []

  for f in (File_Tuple):
    handle = file_handle(f)

    if re.search(r'.sdf', f):
      Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=False)
             if x is not None]

    if re.search(r'.smi', f):
      with open(f, 'r') as fi:
        first_line = fi.readline()
        print(first_line)
     
      if re.search(r'smiles', first_line, re.IGNORECASE):
        Mol = [x for x in Chem.SmilesMolSupplier(handle, titleLine=True,
                 delimiter=' |\t|,') if x is not None]
      else:
        Mol = [x for x in Chem.SmilesMolSupplier(handle, titleLine=False,
                 delimiter=' |\t|,') if x is not None]

    print( "# Found mol in {0}: {1}".format(f, len(Mol)))
    for mol in Mol: List.append(mol)

  gc.collect()
  return List



mols = rdkit_open([sys.argv[1]])

w = Chem.SmilesWriter(sys.argv[2])
for m in mols: w.write(m)

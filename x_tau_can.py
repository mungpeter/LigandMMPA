#!/usr/bin/env python3

import sys
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from molvs import tautomer
#from molvs import standardize
from pathos import multiprocessing


chunk = 500

s_df = pd.read_csv(sys.argv[1], sep='\t', header=0, comment='#')
print(s_df.smiles[:10])

#########################################
class CanonSmiles(object):
  def __init__(self, canon=None):
    self.canon = tautomer.TautomerCanonicalizer()

  def __call__(self, row):
    return self.canonize(row)

  def canonize(self, row ):
    t1 = Chem.MolFromSmiles(row['smiles'])
    r1 = Chem.MolToSmiles(self.canon.canonicalize(t1))
    return r1

##########################################

can = CanonSmiles()
mpi = multiprocessing.Pool()

#s_df['Canon'] = mpi.map(can, s_df, chunk)
s_df['Canon'] = [can(dict(x)) for i, x in tqdm(s_df.iterrows(), total=len(s_df))]
print(len(s_df))

x_df = s_df.drop_duplicates('Canon')
print(len(x_df))
x_df.to_csv(sys.argv[2], sep='\t', header=True, index=False)

########################################################################
#
# Peter M.U. Ung @ Yale/MSSM
#
# v1.0  19.05.16
#
#
# 
#

#!/usr/bin/env python3

import sys,re,os
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole

xxx = Chem.MolFromSmiles('C/C=C\C') ; xxx
xxx


example = r'C[*][*][*]C=C.C[*]CC(=O)OC.c1[*]ccco1.[*]/C(O)=N\C1CC1'
example = r'C[*][*][*]C=C.[*]/C=C\C(=O)OC.c1[*]ccco1.[*]/C(O)=N\C1CC1'
frag = re.sub(r'\\', '{0}', example); frag
Chem.MolFromSmiles(example)

heads = re.findall(r'\[\*\](\/.?)', frag, re.DOTALL) ; heads
for head in heads:
  frag = re.sub(r'\[\*\](\/.?)', head+'[*]', frag)
frag
s = "\\"
s
'{0:s}'.format(s)
print('{:s}'.format(s))



x2 = frag.replace('{0}', '{0}'.format(s))
x2 = frag.format(s)
x2
Chem.MolFromSmiles(x2)
x2 = (frag.format(s)) ; x2
print(x2)

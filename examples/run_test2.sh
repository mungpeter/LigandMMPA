
## parse fragment json into smiles format
../2_parse_mmpdb_frag.py    \
  -list test2_frag.list     \
  -size 25                  \
  -out test.smi \
  -regex 'Se|Te|B'


## number of unique frag:  1641
#  1_abl1.zinc12.cons_10-5.rdkit.fragments.smi.txt  smiles in csv format for step3
#  1_abl1.zinc12.cons_10-5.rdkit.fragments.smi.pkl  smiles in pickle, for step3

#  1_abl1.zinc12.cons_10-5.rdkit.fragments.smi.smi  smiles format, for viewing only

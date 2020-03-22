
## convert to rdkit-compatible smiles format
../0_canonical_smiles_convert.py  \
  1_abl1.zinc12.cons_10-5.smi \
  1_abl1.zinc12.cons_10-5.rdkit.smi

## zip it and put the name in a list
ls 1_abl1.zinc12.cons_10-5.rdkit.smi > step1_smi.list

## generate mmpDB fragments 
../1_mmpdb_frag_gen.csh \
  step1_smi.list        \
  /home/pmung/Dropbox/9_scripts/1_Docking/3_decorate_gen


## result
#  fragment record: 2323
#  1_abl1.zinc12.cons_10-5.fragments.gz

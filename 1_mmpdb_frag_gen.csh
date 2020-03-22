#!/bin/csh

if ($#argv != 2) then
  echo ''
  echo '  > x.csh '
  echo '      [ list of SMILES files ]'
  echo '      [ path to mmpdb-master ]'
  echo ''
endif

echo `pwd` 

foreach smi (`cat $argv[1]`)

  echo $smi
  set name = `basename $smi .smi`

  python $argv[2]/mmpdb-master/mmpdb fragment $smi  \
  -o $name.fragments.json.gz                        \
  -j 10                                             \
  --out fragments.json.gz                           \
  --has-header --max-heavies 67

end


######################################################
#
#   Peter M.U. Ung @ MSSM/Yale
#
#   19.05.10  v1
#
#   Parse chemical library (SMILES) with mmPDB to generate fragments in 
#   JSON format (formatted list of lists)
#
#

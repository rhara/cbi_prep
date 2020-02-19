#!/bin/bash

pdbid=$1
ligname=$2
ncpu=4
npose=4

python extract_lig.py ${pdbid}.pdb.gz ${ligname}
python make_apo.py ${pdbid}.pdb.gz ${ligname}.pdb ${pdbid}_apo.pdb
python write_tleaprc.py ${pdbid} ${pdbid}_apo.pdb ${pdbid}_apo_H.pdb ${pdbid}_apo_H_ref.mol2 > ${pdbid}_tleaprc
tleap -s -f ${pdbid}_tleaprc
obabel ${pdbid}_apo_H.pdb -O ${pdbid}_apo_H_nocharge.mol2 2> /dev/null
python make_charged_protein.py ${pdbid}_apo_H_nocharge.mol2 ${pdbid}_apo_H_ref.mol2 ${pdbid}_apo_H_charged.mol2
smina -r ${pdbid}_apo_H_charged.mol2 -l ${ligname}.pdb --cpu ${ncpu} $(python center.py ${ligname}.pdb) --num_modes ${npose} -o ${ligname}_redock.sdf --log ${pdbid}_smina.log
python rmsd.py ${ligname}.pdb ${ligname}_redock.sdf

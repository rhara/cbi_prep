#!/bin/bash

pdbid=$1
ligname=$2
ncpu=12
npose=4

run () {
    echo $* >> LOG
    $*
}

echo '----------------------' >> LOG
run python extract_lig.py ${pdbid}.pdb.gz ${ligname}
run obabel ${ligname}.pdb -O ${ligname}.sdf
run python make_apo.py ${pdbid}.pdb.gz ${ligname}.pdb ${pdbid}_apo.pdb
run python write_tleaprc.py ${pdbid} ${pdbid}_apo.pdb ${pdbid}_apo_H.pdb ${pdbid}_apo_H_ref.mol2 > ${pdbid}_tleaprc
run tleap -s -f ${pdbid}_tleaprc

run obabel ${pdbid}_apo_H.pdb -O ${pdbid}_apo_H_nocharge.mol2 2> /dev/null
run python make_charged_protein.py ${pdbid}_apo_H_nocharge.mol2 ${pdbid}_apo_H_ref.mol2 ${pdbid}_apo_H_charged.mol2
run smina -r ${pdbid}_apo_H_charged.mol2 -l ${ligname}.sdf --cpu ${ncpu} $(python center.py ${ligname}.pdb) --num_modes ${npose} -o ${ligname}_redock.sdf --seed 0 --log ${pdbid}_smina.log
run python rmsd.py ${ligname}.sdf ${ligname}_redock.sdf

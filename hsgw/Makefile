.SHELL = bash

default:
	@echo Please specify target
	@for f in $$(cat Makefile | grep -E '^\w.*:' | cut -d: -f1) ; do echo '-' $$f ; done

PDB=1mwt
LIG=PNM
test:
	python extract_lig.py $(PDB).pdb.gz $(LIG)
	python make_apo.py $(PDB).pdb.gz $(LIG).pdb $(PDB)_apo.pdb
	python write_tleaprc.py $(PDB) $(PDB)_apo.pdb $(PDB)_apo_H.pdb $(PDB)_apo_H_ref.mol2 > $(PDB)_tleaprc
	tleap -s -f $(PDB)_tleaprc
	obabel $(PDB)_apo_H.pdb -O $(PDB)_apo_H_nocharge.mol2 2> /dev/null
	python make_charged_protein.py $(PDB)_apo_H_nocharge.mol2 $(PDB)_apo_H_ref.mol2 $(PDB)_apo_H_charged.mol2
	smina -r $(PDB)_apo_H_charged.mol2 -l $(LIG).pdb --cpu 4 $$(python center.py $(LIG).pdb) --num_modes 4 -o $(PDB)_redock.sdf --log $(PDB)_smina.log

clean:
	rm -f *~ *.bak *.pdbqt *.log *.pdb *.sdf *.mol2 *_tleaprc
	rm -f *_tleaplog *_obabellog *_rmsd

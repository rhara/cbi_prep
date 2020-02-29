#!/usr/bin/env python

import argparse
import prody
from cbilib import prody_util, ob_util, rdkit_util

def main(args):
    pdbid = args.pdbid
    ligname = args.ligname
    ligand_thres = args.ligand_thres
    pocket_thres = args.pocket_thres

    prody.fetchPDB(pdbid)
    iname = f'{pdbid}.pdb.gz'
    protein = prody.parsePDB(iname, model=1)
    apo = protein.select('protein and not water and not hydrogen')

    ligand = protein.select(f'hetero and resname {ligname} and not hydrogen')
    ligand = prody_util.ligand_pick_one(ligand)
    ligand = ligand.toAtomGroup()
    natoms = ligand.numAtoms()
    print('ligand.numAtoms =', natoms)

    apo = prody_util.apo_near_ligand(apo, ligand, thres=ligand_thres)
    apo = apo.toAtomGroup()

    pocket = prody_util.make_pocket(apo, ligand, thres=pocket_thres)
    pocket = pocket.toAtomGroup()
    nres = pocket.numResidues()
    print('pocket.numResidue =', nres)

    prody.writePDB(f'{pdbid}.apo.pdb.gz', apo)
    prody.writePDB(f'{pdbid}_{ligname}.pdb', ligand)
    prody.writePDB(f'{pdbid}.pocket_{pocket_thres}.pdb', pocket)

    smiles = ob_util.pdb_to_smistring(f'{pdbid}_{ligname}.pdb')
    print('Ligand smiles =', smiles)
    ob_util.pdb_to_sdfile(f'{pdbid}_{ligname}.pdb', f'{pdbid}_{ligname}.sdf', title=f'{pdbid}_{ligname}')
    rdkit_util.make_pair(f'{pdbid}_{ligname}.sdf', f'{pdbid}.pocket_{pocket_thres}.pdb', f'{pdbid}.pair_{pocket_thres}.pkl')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdbid', type=str)
    parser.add_argument('ligname', type=str)
    parser.add_argument('--ligand-thres', type=float, default=5.0)
    parser.add_argument('--pocket-thres', type=float, default=5.0)
    args = parser.parse_args()
    main(args)

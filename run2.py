#!/usr/bin/env python

from openeye.oechem import *
import os, sys, argparse
sys.path.insert(0, os.getcwd())
from cbilib.io import read_pdb, write_molecule
from cbilib.molop import get_ligand_and_protein, make_pocket, split_ligands
import multiprocessing as mp

NOT_RESOLVED = 0
OK = 1
SKIPPED = 2

def worker(args):
    count, pdb_dir, odir, pdbid, ref_ligandname = args
    pro_name = f'{odir}/{pdbid}/{pdbid}_apo.pdb.gz'
    if os.path.exists(pro_name):
        return count, pdbid, SKIPPED, None, None, None, None
    pdbmol = read_pdb(f'{pdb_dir}/{pdbid}.pdb.gz')
    success, ligand, protein = get_ligand_and_protein(pdbmol, ref_ligandname)
    ligands = split_ligands(ligand)
    smis = set()
    for i in range(len(ligands)):
        submol = ligands[i]
        smis.add(OECreateSmiString(submol))
    if 1 < len(smis):
        return count, pdbid, NOT_RESOLVED, ref_ligandname, None, None, None
    ligand = ligands[0]
    pocket = make_pocket(ligand, protein)
    return count, pdbid, OK, ref_ligandname, ligand, protein, pocket

def main(args):
    pdb_dir = args.pdb
    text = 'run1.log'
    odir = args.out

    def gen():
        count = 0
        for line in open(text, 'rt'):
            if 'multiple' in line:
                count += 1
                it = line.split()
                pdbid = it[0]
                ref_ligandname = it[1]
                yield count, pdb_dir, odir, pdbid, ref_ligandname

    log = open('run2.log', 'wt')
    pool = mp.Pool(mp.cpu_count())
    for count, pdbid, ret, ref_ligandname, ligand, protein, pocket in pool.imap_unordered(worker, gen()):
        if ret == OK:
            print(count, pdbid, 'write')
            sdf_name = f'{odir}/{pdbid}/{ref_ligandname}.sdf'
            pro_name = f'{odir}/{pdbid}/{pdbid}_apo.pdb.gz'
            pocket_name = f'{odir}/{pdbid}/{pdbid}_pocket.pdb'
            os.makedirs(f'{odir}/{pdbid}', exist_ok=True)
            write_molecule(ligand, sdf_name)
            write_molecule(protein, pro_name)
            write_molecule(pocket, pocket_name)
        elif ret == SKIPPED:
            print(count, pdbid, 'skipped')
        elif ret == NOT_RESOLVED:
            print(count, pdbid, ref_ligandname, 'error: multiple ligand not resolved')
            print(pdbid, ref_ligandname, 'error: multiple ligand not resolved', file=log)
        else:
            print(count, pdbid, f'error: unknown error ({ret})')
            print(pdbid, f'error: unknown error ({ret})', file=log)
    log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', '-p', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()
    main(args)

#!/usr/bin/env python

from openeye.oechem import *
import argparse, os, sys
sys.path.insert(0, os.getcwd())
from cbilib.io import read_pdb, write_molecule
from cbilib.catalog import get_catalog, get_ligandname_from_catalog
from cbilib.molop import get_ligand_and_protein, make_pocket
import multiprocessing as mp

NOT_FOUND = 0
FOUND = 1
MULTIPLE_LIGAND = 2
SKIPPED = 3

def worker(args):
    count, pdb_dir, odir, pdbid, ref_ligandname = args
    pro_name = f'{odir}/{pdbid}/{pdbid}_apo.pdb.gz'
    if os.path.exists(pro_name):
        return count, pdbid, SKIPPED, None, None, None, None
    pdb_fname = f'{pdb_dir}/{pdbid}.pdb.gz'
    pdbmol = read_pdb(pdb_fname)
    found, ligand, protein = get_ligand_and_protein(pdbmol, ref_ligandname)

    if found:
        ret = FOUND
    else:
        return count, pdbid, NOT_FOUND, ref_ligandname, None, None, None

    smiles = OEMolToSmiles(ligand)

    if '.' in smiles:
        return count, pdbid, MULTIPLE_LIGAND, ref_ligandname, None, None, None

    pocket = make_pocket(ligand, protein)
    return count, pdbid, ret, ref_ligandname, ligand, protein, pocket

def main(args):
    pdb_dir = args.pdb
    csv_name = args.csv
    odir = args.out
    df = get_catalog(csv_name)

    def gen():
        count = 0
        for i in df.index:
            count += 1
            pdbid = df.loc[i, 'PDB_code']
            ref_ligandname = get_ligandname_from_catalog(df, pdbid)
            yield count, pdb_dir, odir, pdbid, ref_ligandname

    log = open('run1.log', 'wt')
    pool = mp.Pool(mp.cpu_count())
    for count, pdbid, ret, ref_ligandname, ligand, protein, pocket in pool.imap_unordered(worker, gen()):
        if ret == FOUND:
            print(count, pdbid, 'write')
            sdf_name = f'{odir}/{pdbid}/{ref_ligandname}.sdf'
            pro_name = f'{odir}/{pdbid}/{pdbid}_apo.pdb.gz'
            pocket_name = f'{odir}/{pdbid}/{pdbid}_pocket.pdb'
            os.makedirs(f'{odir}/{pdbid}', exist_ok=True)
            write_molecule(ligand, sdf_name)
            write_molecule(protein, pro_name)
            write_molecule(pocket, pocket_name)
        elif ret == NOT_FOUND:
            print(count, pdbid, ref_ligandname, 'error: ligand not found')
            print(pdbid, ref_ligandname, 'error: ligand not found', file=log)
        elif ret == SKIPPED:
            print(count, pdbid, 'skipped')
        elif ret == MULTIPLE_LIGAND:
            print(count, pdbid, ref_ligandname, 'error: multiple ligand was found')
            print(pdbid, ref_ligandname, 'error: multiple ligand was found', file=log)
        else:
            print(count, pdbid, f'error: unknown error ({ret})')
            print(pdbid, f'error: unknown error ({ret})', file=log)
    log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', '-p', type=str, required=True)
    parser.add_argument('--csv', '-i', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()
    main(args)

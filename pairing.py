#!/usr/bin/env python

import argparse, os, sys, pickle, gzip
sys.path.insert(0, os.getcwd())
from cbilib.rdkit_io import read_sdf, read_pdb
import multiprocessing as mp

def worker(args):
    pdbid, sdf_name, pocket_name, oname = args
    ligand = read_sdf(sdf_name)
    pocket = read_pdb(pocket_name)
    pickle.dump((ligand, pocket), gzip.open(oname, 'wb'), protocol=4)
    return pdbid, oname

def main(args):
    def gen():
        for root, dirs, files in os.walk(args.dir):
            pdbid = root[-4:]
            sdf_name = None
            pocket_name = None
            oname = root + '/' + pdbid + '_pair.pkl.gz'
            if os.path.exists(oname):
                print(pdbid, 'skipped')
                continue
            for f in files:
                if f.endswith('.sdf'):
                    sdf_name = root + '/' + f
                elif f.endswith('_pocket.pdb'):
                    pocket_name = root + '/' + f
            if sdf_name and pocket_name:
                yield pdbid, sdf_name, pocket_name, oname

    pool = mp.Pool(mp.cpu_count())
    for pdbid, oname in pool.imap_unordered(worker, gen()):
        print(pdbid, '=>', oname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', type=str, required=True)
    args = parser.parse_args()
    main(args)



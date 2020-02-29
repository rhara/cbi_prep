#!/usr/bin/env python

import os, gzip, uuid, argparse, shutil
import openbabel as ob
import subprocess as sp
import multiprocessing as mp
from cbilib import catalog

def prep_protein(uid, protein_iname):
    if protein_iname.endswith('.gz'):
        openf = gzip.open
    else:
        openf = open
    protein_cont = openf(protein_iname, 'rt').read()
    protein_oname = f'{uid}/protein.pdb'
    open(protein_oname, 'wt').write(protein_cont)

def prep_ligand(uid, ligand_iname):
    ligand_cont = open(ligand_iname, 'rt').read()
    ligand_oname = f'{uid}/ligand.sdf'
    open(ligand_oname, 'wt').write(ligand_cont)

def tleap(uid):
    global THISDIR
    cont = open(THISDIR + '/tleaprc.template', 'rt').read()
    cont = cont.replace('{uid}', uid)
    open(f'{uid}/tleaprc', 'wt').write(cont)
    sp.call(f'tleap -s -f {uid}/tleaprc', shell=True)
    sp.call(f'obabel {uid}/protein.mol2 -O {uid}/protein_charged.mol2 > /dev/null 2>&1', shell=True)

def smina(uid, ncpu=None, num_modes=4, seed=0):
    global THISDIR
    if ncpu is None:
        ncpu = os.cpu_count()
    receptor = f'{uid}/protein_charged.mol2'
    ligand = f'{uid}/ligand.sdf'
    oname = f'{uid}/docked.sdf'
    logname = f'{uid}/smina.log'
    boxpars = sp.check_output(f'python {THISDIR}/center.py {ligand}', shell=True).decode().strip()
    sp.call(f'smina -r {receptor} -l {ligand} {boxpars} --cpu {ncpu} --num_modes {num_modes} --seed {seed} -o {oname} --log {logname}', shell=True)

def rmsd(uid):
    global THISDIR
    ref = f'{uid}/ligand.sdf'
    fit = f'{uid}/docked.sdf'
    oname = f'{uid}/rmsd'
    cont = sp.check_output(f'python {THISDIR}/rmsd.py {ref} {fit}', shell=True).decode().strip()
    print(cont)
    open(oname, 'wt').write(cont + '\n')

def retrieve(uid, pdbid):
    global DATADIR
    sp.call(f'cp -av {uid}/smina.log {DATADIR}/{pdbid}/smina.log', shell=True)
    sp.call(f'cp -av {uid}/docked.sdf {DATADIR}/{pdbid}/docked.sdf', shell=True)
    sp.call(f'cp -av {uid}/tleap.log {DATADIR}/{pdbid}/tleap.log', shell=True)
    sp.call(f'cp -av {uid}/rmsd {DATADIR}/{pdbid}/rmsd', shell=True)

THISDIR = None
DATADIR = None

def worker(args):
    global THISDIR

    protein_iname, ligand_iname = args

    print('>', protein_iname)
    print('>', ligand_iname)
    uid = uuid.uuid4().hex
    os.makedirs(uid, mode=0o775)
    prep_protein(uid, protein_iname)
    prep_ligand(uid, ligand_iname)
    tleap(uid)
    smina(uid)
    rmsd(uid)
    pdbid = os.path.basename(protein_iname)[:4]
    retrieve(uid, pdbid)
    shutil.rmtree(uid)

def main(args):
    global THISDIR, DATADIR

    csv_iname = args.csv_iname
    idir = args.idir

    cat = catalog.Catalog(csv_iname)
    print(cat.df)

    THISDIR = os.path.abspath(os.path.dirname(__file__))
    DATADIR = idir

    count = 0
    for i in cat.df.index:
        r = cat.df.loc[i]
        pdbid = r['pdbid']
        ligname = r['ligname']
        protein_iname = f'{idir}/{pdbid}/{pdbid}.apo.pdb.gz'
        ligand_iname = f'{idir}/{pdbid}/{pdbid}_{ligname}.sdf'
        if not (os.path.exists(protein_iname) and os.path.exists(ligand_iname)):
            continue
        count += 1
        print('#'*100)
        print(count, protein_iname, ligand_iname)
        print('#'*100)
        args = (protein_iname, ligand_iname)
        worker(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_iname', type=str)
    parser.add_argument('idir', type=str)
    args = parser.parse_args()
    main(args)

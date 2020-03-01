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
    sp.call(f'tleap -s -f {uid}/tleaprc > /dev/null 2>&1', shell=True)
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
    sp.call(f'smina -r {receptor} -l {ligand} {boxpars} --cpu {ncpu} --num_modes {num_modes} --seed {seed} -o {oname} --log {logname} > /dev/null', shell=True)

def rmsd(uid):
    global THISDIR
    ref = f'{uid}/ligand.sdf'
    fit = f'{uid}/docked.sdf'
    oname = f'{uid}/rmsd'
    cont = sp.check_output(f'python {THISDIR}/rmsd.py {ref} {fit}', shell=True).decode().strip()
    rmsds = []
    for line in cont.split('\n'):
        it = line.split()
        rmsds.append((int(it[0]), round(float(it[1]), 3)))
    open(oname, 'wt').write(cont + '\n')
    return rmsds

def retrieve(uid, pdbid):
    global DATADIR
    sp.call(f'cp -a {uid}/smina.log {DATADIR}/{pdbid}/smina.log', shell=True)
    sp.call(f'cp -av {uid}/docked.sdf {DATADIR}/{pdbid}/docked.sdf', shell=True)
    sp.call(f'cp -a {uid}/tleap.log {DATADIR}/{pdbid}/tleap.log', shell=True)
    sp.call(f'cp -av {uid}/rmsd {DATADIR}/{pdbid}/rmsd', shell=True)

THISDIR = None
DATADIR = None

def worker(args):
    global THISDIR

    count, protein_iname, ligand_iname = args

    print(count, protein_iname, ligand_iname, '[start]')

    pdbid = os.path.basename(protein_iname)[:4]

    uid = uuid.uuid4().hex
    os.makedirs(uid, mode=0o775, exist_ok=True)
    try:
        prep_protein(uid, protein_iname)
        prep_ligand(uid, ligand_iname)
        tleap(uid)
        smina(uid, ncpu=os.cpu_count())
        rmsds = rmsd(uid)
        pdbid = os.path.basename(protein_iname)[:4]
        retrieve(uid, pdbid)
        success = True
    except:
        success = False
        rmsds = None
    shutil.rmtree(uid)
    print(success, count, pdbid, '[done]')
    if rmsds:
        for i, v in rmsds:
            print(count, pdbid, '[rmsd]', i, v)
    return success, count

def main(args):
    global THISDIR, DATADIR

    csv_iname = args.csv_iname
    idir = args.idir
    if args.label:
        q_labels = args.label.split(':')
    else:
        q_labels = None
    cat = catalog.Catalog(csv_iname)

    THISDIR = os.path.abspath(os.path.dirname(__file__))
    DATADIR = idir

    def gen():
        count = 0
        for i in cat.df.index:
            r = cat.df.loc[i]
            pdbid = r['pdbid']
            ligname = r['ligname']
            labels = r['label'].split(':')
            if q_labels:
                ok = True
                for q in q_labels:
                    if q not in labels:
                        ok = False
                        break
                if not ok:
                    continue
            protein_iname = f'{idir}/{pdbid}/{pdbid}.apo.pdb.gz'
            ligand_iname = f'{idir}/{pdbid}/{pdbid}_{ligname}.sdf'
            if not (os.path.exists(protein_iname) and os.path.exists(ligand_iname)):
                continue
            count += 1
            if os.path.exists(f'{idir}/{pdbid}/rmsd'):
                continue
            args = (count, protein_iname, ligand_iname)
            yield args

    pool = mp.Pool(4)
    for ret in pool.imap_unordered(worker, gen()):
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_iname', type=str)
    parser.add_argument('idir', type=str)
    parser.add_argument('--label', '-l', type=str)
    args = parser.parse_args()
    main(args)

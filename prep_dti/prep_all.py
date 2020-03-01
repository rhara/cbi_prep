#!/usr/bin/env python

import os, sys, argparse
import prody
from cbilib import catalog, prody_util, ob_util, rdkit_util
import multiprocessing as mp

csv_name = None
data_dir = None
ligand_thres = None
pocket_thres = None

def test(obj, obj_name, ret):
    r = ret[obj_name + '.success'] = obj is not None
    return r

def worker(args):
    global data_dir, pocket_thres, ligand_thres
    count, pdbid, ligname = args
    os.makedirs(f'{data_dir}/{pdbid}', exist_ok=True, mode=0o755)
    prody_util.fetch_pdb(pdbid, f'{data_dir}/{pdbid}')
    ret = dict(count=count, pdbid=pdbid, ligname=ligname)
    iname = f'{data_dir}/{pdbid}/{pdbid}.pdb.gz'
    protein = prody.parsePDB(iname, model=1)

    ### Temporarily make apo of only protein component
    apo = protein.select('protein and not water and not hydrogen')
    if not test(apo, 'apo', ret): return ret

    ### Temporarily make ligand picked by resname
    ligand = protein.select(f'hetero and resname {ligname} and not hydrogen')
    if not test(ligand, 'lig', ret): return ret

    ### Multiple component will be parsed for one to pick the biggest
    ligand = prody_util.ligand_pick_one(ligand)
    if not test(ligand, 'lig', ret): return ret
    ligand = ligand.toAtomGroup()
    natoms = ligand.numAtoms()
    ret['natoms'] = natoms

    ### Chains limited to near ligand
    apo = prody_util.apo_near_ligand(apo, ligand, thres=ligand_thres)
    if not test(apo, 'apo', ret): return ret
    apo = apo.toAtomGroup()

    ### Residues limited to near ligand
    pocket = prody_util.make_pocket(apo, ligand, thres=pocket_thres)

    if not test(pocket, 'pocket', ret): return ret
    pocket = pocket.toAtomGroup()
    nres = pocket.numResidues()
    ret['nres'] = nres

    prody.writePDB(f'{data_dir}/{pdbid}/{pdbid}.apo.pdb.gz', apo)
    prody.writePDB(f'{data_dir}/{pdbid}/{pdbid}_{ligname}.pdb', ligand)
    prody.writePDB(f'{data_dir}/{pdbid}/{pdbid}.pocket_{pocket_thres}.pdb', pocket)

    smiles = ob_util.pdb_to_smistring(f'{data_dir}/{pdbid}/{pdbid}_{ligname}.pdb')
    ret['smiles'] = smiles
    ob_util.pdb_to_sdfile(f'{data_dir}/{pdbid}/{pdbid}_{ligname}.pdb', f'{data_dir}/{pdbid}/{pdbid}_{ligname}.sdf', title=f'{pdbid}_{ligname}')
    rdkit_util.make_pair(f'{data_dir}/{pdbid}/{pdbid}_{ligname}.sdf', f'{data_dir}/{pdbid}/{pdbid}.pocket_{pocket_thres}.pdb', f'{data_dir}/{pdbid}/{pdbid}.pair_{pocket_thres}.pkl')

    return ret

def main(args):
    global csv_name, data_dir, ligand_thres, pocket_thres
    csv_name = args.csv_name
    data_dir = args.data_dir
    ligand_thres = args.ligand_thres
    pocket_thres = args.pocket_thres
    cat = catalog.Catalog(csv_name)

    os.makedirs(data_dir, exist_ok=True, mode=0o755)

    def gen():
        count = 0
        for i in cat.df.index:
            count += 1
            pdbid = cat.df.loc[i, 'pdbid']
            ligname = cat.get_ligname(pdbid)
            yield count, pdbid, ligname

    pool = mp.Pool(mp.cpu_count())
    for ret in pool.imap_unordered(worker, gen()):
        print(ret['count'], ret['pdbid'], ret['ligname'], ret.get('natoms',None), ret.get('nres',None), end=' ')
        s = [k[:-8] for k in ret if k.endswith('.success') and ret[k]]
        print(','.join(s), end=' ')
        print(ret.get('smiles'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_name', type=str)
    parser.add_argument('data_dir', type=str)
    parser.add_argument('--ligand-thres', '-l', type=float, default=5.0)
    parser.add_argument('--pocket-thres', '-p', type=float, default=5.0)
    args = parser.parse_args()
    main(args)

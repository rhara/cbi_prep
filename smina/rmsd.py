#!/usr/bin/env python

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import math
from rdkit.Chem import MolStandardize


ATOMS = dict(H=1, B=5, C=6, N=7, O=8, F=9, SI=14, P=15, S=16, CL=17, SE=34, BR=35, I=53)

"""
normalizer = MolStandardize.normalize.Normalizer()
lfc = MolStandardize.fragment.LargestFragmentChooser()

def read_pdb_new(fname):
    mol = Chem.MolFromPDBFile(fname, sanitize=False, removeHs=False)
    mol = normalizer.normalize(mol)
    mol = lfc.choose(mol)
    mol = Chem.RemoveHs(mol)
    return mol

def read_sdf_new(fname):
    for mol in Chem.SDMolSupplier(fname, sanitize=False, removeHs=False):
        mol = normalizer.normalize(mol)
        mol = lfc.choose(mol)
        mol = Chem.RemoveHs(mol)
        yield mol
"""

def mk_mol_atom_only(atoms):
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, (el, x, y, z) in enumerate(atoms):
        atom = Chem.Atom(el)
        mol.AddAtom(atom)
        conf.SetAtomPosition(i, [x, y, z])
    mol.AddConformer(conf)
    return Chem.Mol(mol)

def read_sdf(fname):
    f = open(fname, 'rt')
    try:
        while True:
            next(f)
            next(f)
            next(f)
            countline = next(f)
            it = countline.split()
            natoms = int(it[0])
            atoms = []
            for i in range(natoms):
                it = next(f).strip().split()
                x, y, z, el = float(it[0]), float(it[1]), float(it[2]), it[3].upper()
                if el not in ATOMS:
                    raise Exception(f'{el} not in ATOMS')
                if el == 'H':
                    continue
                atoms.append((ATOMS[el], x, y, z))
            yield mk_mol_atom_only(atoms)
            while True:
                line = next(f).strip()
                if line == '$$$$':
                    break
    except StopIteration:
        pass

def get_rmsd(mol1, mol2):
    # if Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2):
    #     return -1.0
    order1 = [x[1] for x in sorted(list(zip(Chem.CanonicalRankAtoms(mol1), range(mol1.GetNumAtoms()))))]
    order2 = [x[1] for x in sorted(list(zip(Chem.CanonicalRankAtoms(mol2), range(mol2.GetNumAtoms()))))]
    conf1 = mol1.GetConformer(0)
    conf2 = mol2.GetConformer(0)
    s = 0
    for i, j in zip(order1, order2):
        c1 = np.array(conf1.GetAtomPosition(i))
        c2 = np.array(conf2.GetAtomPosition(j))
        d2 = np.sum((c1 - c2)*(c1 - c2))
        s += d2
    rmsd = math.sqrt(s/mol1.GetNumAtoms())
    return rmsd

def presanitize(mol):
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        deg = atom.GetDegree()
        if n == 6:
            if 4 < deg:
                atom.SetFormalCharge(deg-4)
        elif n == 7:
            if 3 < deg:
                atom.SetFormalCharge(deg-3)
        elif n == 8:
            if 2 < deg:
                atom.SetFormalCharge(deg-2)

def connect_nearest(mol):
    rwmol = Chem.RWMol(mol)
    conf = rwmol.GetConformer(0)
    coords = np.array(([list(conf.GetAtomPosition(i)) for i in range(rwmol.GetNumAtoms())]))
    thres = 1.9
    thres2 = thres*thres
    connected = set()
    for i in range(coords.shape[0]):
        r = coords[i]
        for j in range(coords.shape[0]):
            if j == i:
                continue
            c = coords[j]
            v = np.sum((c - r)*(c - r))
            if v < thres2:
                if (i, j) in connected:
                    continue
                rwmol.AddBond(i, j, Chem.BondType.SINGLE)
                connected.add((i, j))
                connected.add((j, i))
    mol = Chem.Mol(rwmol)
    presanitize(mol)
    Chem.SanitizeMol(mol)
    return mol

def main(args):
    ref_name = args.ref_name
    fit_name = args.fit_name
    for mol1 in read_sdf(ref_name):
        break
    mol1 = connect_nearest(mol1)
    count = 0
    for mol2 in read_sdf(fit_name):
        count += 1
        mol2 = connect_nearest(mol2)
        rmsd = get_rmsd(mol1, mol2)
        rmsd = round(rmsd, 3)
        print(count, rmsd)
    """
    *Does not work well*
    mol1 = read_pdb_new(ref_name)
    count = 0
    for mol2 in read_sdf_new(fit_name):
        count += 1
        rmsd = get_rmsd(mol1, mol2)
        rmsd = round(rmsd, 3)
        print(count, rmsd)
    """

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_name', type=str)
    parser.add_argument('fit_name', type=str)
    args = parser.parse_args()
    main(args)


from openeye.oechem import *
import os, argparse, pickle
import numpy as np

def read_molecule(fname):
    ifs = oemolistream(fname)
    pdbmol = OEGraphMol()
    OEReadMolecule(ifs, pdbmol)
    return pdbmol

def write_molecule(mol, oname):
    ofs = oemolostream(oname)
    OEWriteMolecule(ofs, mol)
    ofs.close()

def main(args):
    pdb_name = args.pdb
    charge_name = args.charge
    pdbid = os.path.basename(pdb_name)[:4]

    protein = read_molecule(pdb_name)
    charge = read_molecule(charge_name)

    protein_atoms = list(protein.GetAtoms())
    charge_atoms = list(charge.GetAtoms())

    for atom in protein.GetAtoms():
        i = atom.GetIdx()
        ch = charge.GetAtom(OEHasAtomIdx(i)).GetPartialCharge()
        atom.SetPartialCharge(ch)

    basedir = os.path.dirname(args.pdb)
    oname = f'{basedir}/{pdbid}_apo_H_charged.mol2'
    write_molecule(protein, oname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('charge', type=str)
    args = parser.parse_args()
    main(args)

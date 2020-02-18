import sys, os, argparse
from collections import namedtuple
from rdkit import Chem
import numpy as np

ELEMENTS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53}

def get_lines_from_mol2(fname, no=1):
    Line = namedtuple('Line', ['no', 'section', 'linecount', 'cont'])
    Atom = namedtuple('Atom', ['idx', 'name', 'x', 'y', 'z', 'el', 'hyb', 'gr', 'gr_name', 'charge'])
    Bond = namedtuple('Bond', ['idx', 'atom1', 'atom2', 'degree'])
    lines = []
    no = 0
    count = 0
    for line in open(fname)
        line = line.strip()
        if line.startswith('@<TRIPOS>'):
            section = line[9:]
            count = 0
            if section == 'MOLECULE':
                no += 1
            continue
        count += 1
        lines.append(Line(no, section, count, line))
    atoms = []
    for line in lines:
        if line.no != no:
            continue
        if line.section == 'ATOM':
            it = line.cont.split()
            idx = int(it[0])-1
            name = it[1]
            x = float(it[2])
            y = float(it[3])
            z = float(it[4])
            t = it[5].split('.')
            hyb = None
            if 1 < len(t):
                el, hyb = t
            else:
                el = t[0]
            gr = int(it[6])
            gr_name = it[7]
            charge = float(it[8])
            atoms.append(Atom(idx, name, x, y, z, el, hyb, gr, gr_name, charge))
    bonds = []
    for line in lines:
        if 2 <= line.no:
            continue
        if line.section == 'BOND':
            it = line.cont.split()
            idx = int(it[0]) - 1
            atom1 = int(it[1]) - 1
            atom2 = int(it[2]) - 1
            degree = it[3]
            bonds.append(Bond(idx, atom1, atom2, degree))
    return atoms, bonds

def make_atoms_collection(atoms):
    """
    https://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
    """
    rwmol = Chem.RWMol()
    for atom in atoms:
        a = Chem.Atom(ELEMENTS[atom.el])
        if atom.el == 'C':
            if atom.hyb == '3':
                a.SetHybridization(Chem.HybridizationType.SP3)
            elif atom.hyb in ['2', 'ar']:
                a.SetHybridization(Chem.HybridizationType.SP2)
            elif atom.hyb == '1':
                a.SetHybridization(Chem.HybridizationType.SP)
            elif atom.hyb == 'cat':
                a.SetHybridization(Chem.HybridizationType.SP2)
                a.SetFormalCharge(1)
            else:
                raise Exception((sys.argv, atom))
        elif atom.el == 'N':
            if atom.hyb == '4':
                a.SetHybridization(Chem.HybridizationType.SP3)
                a.SetFormalCharge(1)
            elif atom.hyb in ['3', 'am', 'pl3']:
                a.SetHybridization(Chem.HybridizationType.SP3)
            elif atom.hyb in ['2', 'ar']:
                a.SetHybridization(Chem.HybridizationType.SP2)
            elif atom.hyb == '1':
                a.SetHybridization(Chem.HybridizationType.SP)
            else:
                raise Exception((sys.argv, atom))
        elif atom.el == 'O':
            if atom.hyb == 'co2' and atom.name == 'OXT':
                a.SetHybridization(Chem.HybridizationType.SP3)
                a.SetFormalCharge(-1)
            elif atom.hyb == 'co2':
                a.SetHybridization(Chem.HybridizationType.SP2)
            elif atom.hyb == '2':
                a.SetHybridization(Chem.HybridizationType.SP2)
            elif atom.hyb == '3':
                a.SetHybridization(Chem.HybridizationType.SP3)
                a.SetHybridization(Chem.HybridizationType.SP2)
            else:
                raise Exception((sys.argv, atom))
        elif atom.el == 'P':
            if atom.hyb == '3':
                a.SetHybridization(Chem.HybridizationType.SP3)
            else:
                raise Exception((sys.argv, atom))
        elif atom.el == 'S':
            if atom.hyb == '3':
                a.SetHybridization(Chem.HybridizationType.SP3)
            elif atom.el == 'S' and atom.hyb == '2':
                a.SetHybridization(Chem.HybridizationType.SP2)
            elif atom.el == 'S' and atom.hyb == 'O2':
                pass
            else:
                raise Exception((sys.argv, atom))
        elif atom.el == 'H':
            if atom.hyb is None:
                a.SetHybridization(Chem.HybridizationType.S)
        elif atom.el in ['F', 'Cl', 'Br', 'I']:
            if atom.hyb is None:
                pass
            else:
                raise Exception((sys.argv, atom))
        else:
            raise Exception((sys.argv, atom))
        rwmol.AddAtom(a)
    return rwmol

def make_conformer(rwmol, atoms):
    conf = Chem.Conformer(len(atoms))
    for atom in atoms:
        conf.SetAtomPosition(atom.idx, [atom.x, atom.y, atom.z])
    rwmol.AddConformer(conf)

def make_bonding(rwmol, bonds):
    for bond in bonds:
        degree = bond.degree
        btype = None
        if degree == '1':
            btype = Chem.BondType.SINGLE
        elif degree == '2':
            btype = Chem.BondType.DOUBLE
        elif degree == '3':
            btype = Chem.BondType.TRIPLE
        elif degree == 'am':
            btype = Chem.BondType.SINGLE
        elif degree == 'ar':
            btype = Chem.BondType.AROMATIC
        else:
            raise Exception((sys.argv, bond))
        rwmol.AddBond(bond.atom1, bond.atom2, btype)

def make_rdkit_molecule_from_mol2(fname):
    atoms, bonds = get_lines_from_mol2(sys.argv[1])
    rwmol = make_atoms_collection(atoms)
    make_conformer(rwmol, atoms)
    make_bonding(rwmol, bonds)
    mol = Chem.RemoveHs(rwmol, sanitize=False)
    return mol

def get_coords(mol, confidx=0):
    conf = mol.GetConformer(confidx)
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    return coords

def main(args):
    iname = args.iname
    mol = make_rdkit_molecule_from_mol2(iname)
    smiles = Chem.MolToSmiles(mol)
    print(os.path.basename(iname), smiles)
    #coords = get_coords(mol)
    #print(coords)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iname', type=str)
    args = parser.parse_args()
    main(args)

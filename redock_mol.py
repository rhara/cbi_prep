import sys
from collections import namedtuple
from rdkit import Chem

Line = namedtuple('Line', ['no', 'section', 'linecount', 'cont'])
Atom = namedtuple('Atom', ['idx', 'name', 'x', 'y', 'z', 'el', 'hyb', 'gr', 'gr_name', 'charge'])
Bond = namedtuple('Bond', ['idx', 'atom1', 'atom2', 'degree'])
ELEMENTS = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17, 'Br': 35, 'I': 53}
lines = []
no = 0
count = 0
for line in open(sys.argv[1]):
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
    if 2 <= line.no:
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
#help(Chem.RWMol)

mol = Chem.RWMol()
for atom in atoms:
    print(atom)
    a = Chem.Atom(ELEMENTS[atom.el])
    if atom.el == 'C' and atom.hyb == '3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'S' and atom.hyb == '3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'N' and atom.hyb == '4':
        a.SetHybridization(Chem.HybridizationType.SP3)
        a.SetFormalCharge(1)
    elif atom.el == 'H' and atom.hyb is None:
        a.SetHybridization(Chem.HybridizationType.S)
    elif atom.el == 'C' and atom.hyb == '2':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'O' and atom.hyb == 'co2' and atom.name == 'OXT':
        a.SetHybridization(Chem.HybridizationType.SP3)
        a.SetFormalCharge(-1)
    elif atom.el == 'O' and atom.hyb == 'co2':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'N' and atom.hyb == 'am':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'O' and atom.hyb == '2':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'C' and atom.hyb == 'ar':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'O' and atom.hyb == '3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'N' and atom.hyb == 'pl3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'N' and atom.hyb == '2':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'P' and atom.hyb == '3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'C' and atom.hyb == '1':
        a.SetHybridization(Chem.HybridizationType.SP)
    elif atom.el == 'N' and atom.hyb == '3':
        a.SetHybridization(Chem.HybridizationType.SP3)
    elif atom.el == 'N' and atom.hyb == 'ar':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'F' and atom.hyb is None:
        pass
    elif atom.el == 'C' and atom.hyb == 'cat':
        a.SetHybridization(Chem.HybridizationType.SP2)
        a.SetFormalCharge(1)
    elif atom.el == 'S' and atom.hyb == 'O2':
        pass
    elif atom.el == 'S' and atom.hyb == '2':
        a.SetHybridization(Chem.HybridizationType.SP2)
    elif atom.el == 'Br' and atom.hyb is None:
        pass
    elif atom.el == 'I' and atom.hyb is None:
        pass
    elif atom.el == 'Cl' and atom.hyb is None:
        pass
    elif atom.el == 'N' and atom.hyb == '1':
        a.SetHybridization(Chem.HybridizationType.SP)
    else:
        raise Exception((sys.argv, atom))
    mol.AddAtom(a)

for bond in bonds:
    print(bond)
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
    mol.AddBond(bond.atom1, bond.atom2, btype)

mol = Chem.RemoveHs(mol, sanitize=False)
print(Chem.MolToSmiles(mol))

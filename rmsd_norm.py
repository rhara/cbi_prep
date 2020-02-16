from openeye.oechem import *
import sys

def normalize_mol(mol):
    mol = OEGraphMol(mol)
    OESuppressHydrogens(mol)
    strategy = OEUncolorStrategy_RemoveAtomProperties | OEUncolorStrategy_ConvertBondTypeToSingle
    mol_ = OEGraphMol(mol)
    OEUncolorMol(mol_, mol, strategy)
    return mol_

OEThrow.SetLevel(OEErrorLevel_Fatal)
ifs = oemolistream(sys.argv[1])
refmol = OEGraphMol()
OEReadMolecule(ifs, refmol)
ifs.close()

refmol_ = normalize_mol(refmol)

ifs = oemolistream(sys.argv[2])
count = 0
for mol in ifs.GetOEGraphMols():
    count += 1
    mol_ = normalize_mol(mol)
    rmsd = OERMSD(refmol_, mol_, True, True, False)
    rmsd = round(rmsd, 3)
    print(count, rmsd)


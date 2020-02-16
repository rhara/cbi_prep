from openeye.oechem import *
import sys

ifs = oemolistream(sys.argv[1])
refmol = OEGraphMol()
OEReadMolecule(ifs, refmol)
ifs.close()

ifs = oemolistream(sys.argv[2])
count = 0
for mol in ifs.GetOEGraphMols():
    count += 1
    rmsd = OERMSD(refmol, mol, True, True, False)
    rmsd = round(rmsd, 3)
    print(count, rmsd)
ifs.close()

from rdkit import Chem
from rdkit.Chem import MolStandardize
import argparse, os, gzip, sys
import numpy as np
import pandas as pd

def main(args):
    pdb_name = args.pdb_name
    ligand_name = args.ligand_name
    pdbid = os.path.basename(pdb_name)[:4]
    oname = f'{pdbid}_{ligand_name}.pdb'
    openf = gzip.open if pdb_name.endswith('.gz') else open

    atomnames = set()
    ligands = {}
    found = False
    for line in openf(pdb_name, 'rt'):
        line = line.rstrip()
        try:
            pdb_tag = line[:6]
            pdb_ligname = line[17:20].strip()
            pdb_chainId = line[21]
            pdb_atomname = line[13:16].strip()
        except:
            continue
        if pdb_tag == 'HETATM' and pdb_ligname == ligand_name:
            found = True
            if pdb_atomname in atomnames:
                continue
            if pdb_chainId not in ligands:
                ligands[pdb_chainId] = []
            ligands[pdb_chainId].append(line)
            atomnames.add(pdb_atomname)
    if not found:
        return -1

    chainIds = sorted(ligands)
    chainId = chainIds[0]

    mol = Chem.MolFromPDBBlock('\n'.join(ligands[chainId]))
    if mol is None:
        return -1

    print(Chem.MolToSmiles(mol), ligand_name, file=sys.stderr)

    """
    if multiple molecule is found in the set, take the largest one
    """
    normalizer = MolStandardize.normalize.Normalizer()
    _mol = normalizer.normalize(mol)
    lfc = MolStandardize.fragment.LargestFragmentChooser()
    mol = lfc.choose(_mol)
    # mols = list(Chem.GetMolFrags(mol, asMols=True))
    # if 1 < len(mols):
    #     mols.sort(key=lambda m: m.GetNumAtoms(), reverse=True)
    # mol = mols[0]

    lines = Chem.MolToPDBBlock(mol).split('\n')
    with open(oname, 'wt') as out:
        for line in lines:
            if line.startswith('HETATM'):
                out.write(line + '\n')

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_name', type=str)
    parser.add_argument('ligand_name', type=str)
    args = parser.parse_args()
    sys.exit(main(args))

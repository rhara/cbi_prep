from rdkit import Chem
import argparse, os, gzip, sys
import numpy as np

def main(args):
    try:
        pdb_name = args.pdb_name
        ligand_name = args.ligand_name
        pdbid = os.path.basename(pdb_name)[:4]
        oname = f'{pdbid}_{ligand_name}.pdb'
        openf = gzip.open if pdb_name.endswith('.gz') else open

        atomnames = set()
        ligands = {}
        for line in openf(pdb_name, 'rt'):
            line = line.rstrip()
            try:
                pdb_tag = line[:6]
                pdb_ligname = line[17:20]
                pdb_chainId = line[21]
                pdb_atomname = line[13:16].strip()
            except:
                continue
            if pdb_tag == 'HETATM' and pdb_ligname == ligand_name:
                if pdb_atomname in atomnames:
                    continue
                if pdb_chainId not in ligands:
                    ligands[pdb_chainId] = []
                ligands[pdb_chainId].append(line)
                atomnames.add(pdb_atomname)

        chainIds = sorted(ligands)
        chainId = chainIds[0]
        mol = Chem.MolFromPDBBlock('\n'.join(ligands[chainId]))
        print(Chem.MolToSmiles(mol), ligand_name)
        mols = list(Chem.GetMolFrags(mol, asMols=True))
        if 1 < len(mols):
            mols.sort(key=lambda m: m.GetNumAtoms(), reverse=True)
        mol = mols[0]
        lines = Chem.MolToPDBBlock(mol).split('\n')
        #ref_xyzs = np.array([(float(line[30:38]), float(line[38:46]), float(line[46:54])) for line in lines if line.startswith('HETATM')])
        ##print(ref_xyzs)
        ##print('\n'.join(lines))
        with open(oname, 'wt') as out:
            for line in lines:
                if line.startswith('HETATM'):
                    out.write(line + '\n')
    except:
        print(locals(), file=sys.stderr)
        return 1
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_name', type=str)
    parser.add_argument('ligand_name', type=str)
    args = parser.parse_args()
    sys.exit(main(args))

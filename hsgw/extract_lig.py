import argparse
import gzip

"""
0         1         2         3         4         5         6
0123456789012345678901234567890123456789012345678901234567890123456789
HETATM10176  O8  PNM B 708     -33.466  47.734  66.394  0.97 48.08           O
"""
def main(args):
    pdb_name = args.pdb_name
    ligand_name = args.ligand_name
    oname = f'{ligand_name}.pdb'
    openf = gzip.open if pdb_name.endswith('.gz') else open

    ligands = {}
    for line in openf(pdb_name, 'rt'):
        line = line.rstrip()
        try:
            pdb_tag = line[:6]
            pdb_ligname = line[17:20]
            pdb_chainId = line[21]
        except:
            continue
        if pdb_tag == 'HETATM' and pdb_ligname == ligand_name:
            if pdb_chainId not in ligands:
                ligands[pdb_chainId] = []
            ligands[pdb_chainId].append(line)

    chainIds = sorted(ligands)
    chainId = chainIds[0]
    with open(oname, 'wt') as out:
        for line in ligands[chainId]:
            out.write(line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_name', type=str)
    parser.add_argument('ligand_name', type=str)
    args = parser.parse_args()
    main(args)

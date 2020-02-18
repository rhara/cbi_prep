import os, argparse
from rdkit import Chem

def main(args):
    nocharge_name = args.nocharge
    charge_ref = args.charge_ref
    oname = args.oname

    ok = False
    pch = []
    for line in open(charge_ref, 'rt'):
        line = line.strip()
        if line == '@<TRIPOS>ATOM':
            ok = True
            continue
        if line == '@<TRIPOS>BOND':
            ok = False
            continue
        if ok:
            it = line.split()
            pch.append(float(it[8]))

    out = open(oname, 'wt')
    ok = False
    for line in open(nocharge_name):
        line = line.rstrip()
        if line == '@<TRIPOS>ATOM':
            ok = True
            out.write(line + '\n')
            continue
        if line.startswith('@') and line != '@<TRIPOS>ATOM':
            ok = False
            out.write(line + '\n')
            continue
        if ok:
            v = pch.pop(0)
            s = f'{v:7.4f}'
            out.write(line[:69] + s + line[76:] + '\n')
        else:
            out.write(line + '\n')

    # basedir = os.path.dirname(args.pdb)
    # oname = f'{basedir}/{pdbid}_apo_H_charged.mol2'
    # write_molecule(protein, oname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('nocharge', type=str)
    parser.add_argument('charge_ref', type=str)
    parser.add_argument('oname', type=str)
    args = parser.parse_args()
    main(args)

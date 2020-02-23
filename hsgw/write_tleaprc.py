#!/usr/bin/env python

import argparse

def main(args):
    prefix = args.prefix
    pdb_name = args.pdb
    opdb_name = args.opdb
    omol2_name = args.omol2
    rc = open('tleaprc.template').read()
    rc = rc.replace('{prefix}', prefix).\
            replace('{pdb_name}', pdb_name).\
            replace('{opdb_name}', opdb_name).\
            replace('{omol2_name}', omol2_name)
    print(rc)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix', type=str)
    parser.add_argument('pdb', type=str)
    parser.add_argument('opdb', type=str)
    parser.add_argument('omol2', type=str)
    args = parser.parse_args()
    main(args)

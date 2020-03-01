import pandas as pd
import re, sys, math
from collections import namedtuple


def get_kinase_pdbids():
    df = pd.read_csv('KLIFS_export.csv')
    return set(df['PDB'])

def get_index_data(iname, kinases, src):
    AMINO_ACIDS = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()
    UNIT = {'M': 1, 'mM': 1e-3, 'uM': 1e-6, 'nM': 1e-9, 'pM': 1e-12, 'fM': 1e-15}
    value_pat = re.compile('^(IC50|Ki|Kd)(=|>|>=|<|<=|~)([.0-9]+)(M|mM|uM|nM|pM|fM)$')
    Data = namedtuple('Data', ['pdbid', 'ligname', 'measure', 'pvalue', 'src', 'label'])

    data = {}
    for line in open(iname, 'rt'):
        line = line.rstrip()
        if line.startswith('#') or line == '':
            continue
        it = line.split(maxsplit=7)
        pdbid, year, value, ligname = it[0], int(it[2]), it[4], it[7][1:-1]
        if len(ligname) != 3 or ligname.startswith('_'):
            continue
        if ligname in AMINO_ACIDS:
            continue
        m = value_pat.match(value)
        assert m
        measure = m.group(1)
        value = float(m.group(3))
        unit = m.group(4)
        value *= UNIT[unit]
        pvalue = round(-math.log10(value), 3)
        label = 'kinase' if pdbid in kinases else ''
        data[pdbid] = Data(pdbid, ligname, measure, pvalue, src, label)
    return data

def main():
    kinases = get_kinase_pdbids()
    data_2018 = get_index_data('INDEX_general_PL_data.2018', kinases, '2018')
    data_2019 = get_index_data('INDEX_general_PL_data.2019', kinases, '2019')
    keys_2018 = set(data_2018.keys())
    keys_2019 = set(data_2019.keys())
    data_all = []
    for k in keys_2018.union(keys_2019):
        if k in data_2018 and k in data_2019:
            if data_2018[k].pvalue == data_2019[k].pvalue:
                data_all.append(data_2018[k])
            else:
                print('*'*10, k, 'updated in 2019')
        elif k in data_2018 and k not in data_2019:
            data_all.append(data_2018[k])
        elif k not in data_2018 and k in data_2019:
            data_all.append(data_2019[k])
        else:
            print(k, 'impossible')
    df = pd.DataFrame(data_all)
    df.to_csv('data.csv', index=False)

if __name__ == '__main__':
    main()

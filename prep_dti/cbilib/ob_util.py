import openbabel as ob

pdb_to_smi_conv = ob.OBConversion()
pdb_to_smi_conv.SetInAndOutFormats('pdb', 'smi')

pdb_to_sdf_conv = ob.OBConversion()
pdb_to_sdf_conv.SetInAndOutFormats('pdb', 'sdf')

def pdb_to_smistring(iname):
    mol = ob.OBMol()
    pdb_to_smi_conv.ReadFile(mol, iname)
    smiles = pdb_to_smi_conv.WriteString(mol).strip().split()[0]
    return smiles

def pdb_to_sdfile(iname, oname, title=None):
    mol = ob.OBMol()
    pdb_to_smi_conv.ReadFile(mol, iname)
    cont = pdb_to_sdf_conv.WriteString(mol)
    lines = cont.strip().split('\n')
    if title:
        lines[0] = title
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('M') and line.endswith('END'):
            end_mark = i
            break
    del lines[end_mark+1:]
    lines.append('$$$$')
    cont = '\n'.join(lines) + '\n'
    open(oname, 'wt').write(cont)

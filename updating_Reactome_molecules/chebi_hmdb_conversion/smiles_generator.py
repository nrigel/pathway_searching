import openbabel, pybel, pickle

def remove_proton_charges(SDF, DB):
    ID = {'HMDB': 'DATABASE_ID', 'CHEBI': 'ChEBI ID'}
    # remove proton-related charges for better NMR-centric searching
    smiles = {}
    for chebi in pybel.readfile("sdf", SDF):
        can = chebi.write('can')
        if chebi.OBMol.GetTotalCharge() or '+' in can or '-' in can:
            chebi.OBMol.BeginModify()
            chebi.OBMol.SetTotalCharge(0)
            for atom in openbabel.OBMolAtomIter(chebi.OBMol):
                charge = atom.GetFormalCharge()
                if charge > 0:
                    if atom.ImplicitHydrogenCount():
                        atom.SetFormalCharge(0)
                if charge < 0:
                    atom.SetFormalCharge(0)
            chebi.OBMol.EndModify()
        can = chebi.write('can').split()
        if can:
            smiles[chebi.data[ID[DB]]] = can[0]
    return smiles

with open('HMDB_SMILES_noH.pkl', 'wb') as pkl:
    pickle.dump(remove_proton_charges("HMDB_structures.sdf", 'HMDB'), pkl, pickle.HIGHEST_PROTOCOL)

with open('CHEBI_SMILES_noH.pkl', 'wb') as pkl:
    pickle.dump(remove_proton_charges("ChEBI_complete.sdf", 'CHEBI'), pkl, pickle.HIGHEST_PROTOCOL)

# for fast searching... dont need
def fingerprint_sdf(SDF):
    conv = openbabel.OBConversion()
    conv.OpenInAndOutFiles(SDF, SDF[:-4]+'_index.fs')
    conv.SetInAndOutFormats("sdf","fs")
    conv.Convert()




import csv, sys, openbabel, pybel, pickle

sys.path.append('../../scripts')
from neo4j_connect import NMRTServer

NMRT = NMRTServer()
with NMRT._driver.session() as db:
    compounds = db.run('MATCH (m:metabolites) RETURN m.Compound').value()


db_formatted = {}
for compound in compounds:
    formatted = compound
    
    if '_' in compound:
        formatted = compound.split('_')[0]

    if 'HMDB' in formatted:
        formatted = 'HMDB'+formatted[4:].zfill(7)
    
    db_formatted[compound] = formatted   

with open('../HMDB_SMILES_noH.pkl', 'rb') as pkl:
    hmdb_smiles = pickle.load(pkl)
smiles_hmdb = {}
for hmdb in hmdb_smiles:
    if hmdb_smiles[hmdb] not in smiles_hmdb:
        smiles_hmdb[hmdb_smiles[hmdb]] = []
    smiles_hmdb[hmdb_smiles[hmdb]].extend([hmdb])
    
conv = openbabel.OBConversion()
colmar_to_hmdb, BMRB_SMILES = {}, {}
for compound in db_formatted:
    if 'bmse' in db_formatted[compound]:
        iso_smiles, can_smiles = None, None
        with open('bmrb_entries/'+db_formatted[compound]+'/'+db_formatted[compound]+'.str') as dbinfo:
            for line in dbinfo:
                if 'isomeric' in line:
                    if len(line.rstrip().split()) > 1:
                        iso_smiles = line.rstrip().split()[1]
                if 'canonical' in line:
                    can_smiles = line.rstrip().split()[1]
        
        def find_match(smiles):
            bmrb = pybel.readstring('smi', smiles)
            can = bmrb.write('can')
            if bmrb.OBMol.GetTotalCharge() or '+' in can or '-' in can:
                bmrb.OBMol.BeginModify()
                bmrb.OBMol.SetTotalCharge(0)
                for atom in openbabel.OBMolAtomIter(bmrb.OBMol):
                    charge = atom.GetFormalCharge()
                    if charge > 0:
                        if atom.ImplicitHydrogenCount():
                            atom.SetFormalCharge(0)
                    if charge < 0:
                        atom.SetFormalCharge(0)
                bmrb.OBMol.EndModify()

            can = bmrb.write('can').split()
            if can:
                if can[0] in smiles_hmdb:
                    colmar_to_hmdb[compound] = smiles_hmdb[can[0]]
                else:
                    return can[0]

        if iso_smiles:
            smiles = find_match(iso_smiles)
            if smiles:
                BMRB_SMILES[compound] = smiles
        if can_smiles and not iso_smiles:
            smiles = find_match(can_smiles)
            if smiles:
                BMRB_SMILES[compound] = smiles
            
    else:
        colmar_to_hmdb[compound] = db_formatted[compound]

with open('colmar_to_hmdb.pkl', 'wb') as pkl:
    pickle.dump(colmar_to_hmdb, pkl, pickle.HIGHEST_PROTOCOL)

with open('BMRB_SMILES_noH.pkl', 'wb') as pkl:
    pickle.dump(BMRB_SMILES, pkl, pickle.HIGHEST_PROTOCOL)
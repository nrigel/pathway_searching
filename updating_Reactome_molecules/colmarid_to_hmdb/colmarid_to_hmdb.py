import csv, sys, openbabel, pybel, pickle

sys.path.append('../scripts')
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
    
conv = openbabel.OBConversion()
bmrb_to_hmdb = {}
for compound in db_formatted:
    if 'bmse' in db_formatted[compound]:
        with open('bmrb_entries/'+db_formatted[compound]+'/'+db_formatted[compound]+'.str') as dbinfo:
            for line in dbinfo:
                if 'canonical' in line:
                    smiles = line.rstrip().split()[1]
                    break

        bmrb_to_hmdb[compound] = []
        
        conv.OpenInAndOutFiles("../chebi_hmdb_conversion/HMDB_index.fs", "results.sdf")
        conv.SetInAndOutFormats("fs", "sdf")
        conv.AddOption("s", conv.GENOPTIONS, smiles)
        conv.AddOption("e", conv.INOPTIONS)
        conv.Convert()
        
        for hmdb in pybel.readfile("sdf", "results.sdf"):
            bmrb_to_hmdb[compound].extend([hmdb.data['HMDB_ID']])
    else:
        bmrb_to_hmdb[compound] = db_formatted[compound]

with open('bmrb_to_hmdb.pkl', 'wb') as pkl:
    pickle.dump(bmrb_to_hmdb, pkl, pickle.HIGHEST_PROTOCOL)
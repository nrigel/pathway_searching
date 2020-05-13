import pickle, sys
sys.path.append('../scripts')
from neo4j_connect import ReactomeServer
Reactome = ReactomeServer()

with open('colmar_to_hmdb.pkl', 'rb') as pkl:
    colmar_to_hmdb = pickle.load(pkl)

with open('CHEBI_TO_HMDB.pkl', 'rb') as pkl:
    chebi_to_hmdb = pickle.load(pkl)

# want to add COLMAR IDs and HMDB IDs to Reactome molecules
hmdb_to_colmar = {} # reverse dict
for colmar in colmar_to_hmdb:
    if type(colmar_to_hmdb[colmar]) == str:
        if colmar_to_hmdb[colmar] not in hmdb_to_colmar:
            hmdb_to_colmar[colmar_to_hmdb[colmar]] = []
        hmdb_to_colmar[colmar_to_hmdb[colmar]].extend([colmar])
    if type(colmar_to_hmdb[colmar]) == list:
        for hmdb in colmar_to_hmdb[colmar]:
            if hmdb not in hmdb_to_colmar:
                hmdb_to_colmar[hmdb] = []
            hmdb_to_colmar[hmdb].extend([colmar])

chebi_to_colmar = {}
for chebi in chebi_to_hmdb:
    chebi_to_colmar[chebi] = []
    for hmdb in chebi_to_hmdb[chebi]:
        if hmdb in hmdb_to_colmar:
            chebi_to_colmar[chebi].extend(hmdb_to_colmar[hmdb])   

with Reactome._driver.session() as db:
    reactome_compounds = db.run('MATCH (m:ReferenceMolecule) RETURN m.displayName').value()
    for compound in reactome_compounds:
        chebi = ''
        for i in range(compound.find('ChEBI'), len(compound)):
            if compound[i] == ']' or compound[i] == ' ':
                break
            chebi = chebi+compound[i]
        chebi = 'CHEBI:'+chebi.split(':')[1]
        if chebi in chebi_to_colmar: # means that this molecule matched to at least one HMDB molecule
            hmdb = '['+', '.join(['"'+h+'"' for h in chebi_to_hmdb[chebi]])+']'
            colmar = '['+', '.join(['"'+c+'"' for c in chebi_to_colmar[chebi]])+']'
            db.run('MATCH (m:ReferenceMolecule) WHERE m.displayName = "'+compound+'" SET m.HMDB = '+hmdb+', m.COLMAR = '+colmar)
            
    
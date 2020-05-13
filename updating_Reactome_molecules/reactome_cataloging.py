import sys, openbabel, pybel, pickle

sys.path.append('../scripts')
from neo4j_connect import ReactomeServer
Reactome = ReactomeServer()

sys.path.append('/Users/nickrigel/Documents/GitHub/nmr-transformer/nmrt')
from chemstructure import Compound

reactome_to_chebi = {}
with Reactome._driver.session() as db:
    reactome_compounds = db.run('MATCH (m:ReferenceMolecule) RETURN m.displayName').value()

    for compound in reactome_compounds:
        chebi = ''
        for i in range(compound.find('ChEBI'), len(compound)):
            if compound[i] == ']' or compound[i] == ' ':
                break
            chebi = chebi+compound[i]
        chebi = 'CHEBI:'+chebi.split(':')[1]

        reactome_to_chebi[compound] = chebi

with open('HMDB_SMILES_noH.pkl', 'rb') as pkl:
    hmdb_smiles = pickle.load(pkl)
smiles_hmdb = {}
for hmdb in hmdb_smiles:
    if hmdb_smiles[hmdb] not in smiles_hmdb:
        smiles_hmdb[hmdb_smiles[hmdb]] = []
    smiles_hmdb[hmdb_smiles[hmdb]].extend([hmdb])

with open('CHEBI_SMILES_noH.pkl', 'rb') as pkl:
    chebi_smiles = pickle.load(pkl)

with open('BMRB_SMILES_noH.pkl', 'rb') as pkl:
    bmrb_smiles = pickle.load(pkl)
smiles_bmrb = {}
for bmrb in bmrb_smiles:
    if bmrb_smiles[bmrb] not in smiles_bmrb:
        smiles_bmrb[bmrb_smiles[bmrb]] = []
    smiles_bmrb[bmrb_smiles[bmrb]].extend([bmrb])

with open('colmar_to_hmdb.pkl', 'rb') as pkl:
    colmar_to_hmdb = pickle.load(pkl)
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

with Reactome._driver.session() as db:
    for compound in reactome_to_chebi:
        
        chebi = reactome_to_chebi[compound]

        if chebi not in chebi_smiles:
            continue
        
        #print(compound)

        S = Compound(chebi_smiles[chebi])
        smiles = chebi_smiles[chebi]

        if not S.spinsystems(): # skip molecules that don't have carbons because they are irrelavant to us
            continue
        
        S.canonicalize()
        can_smiles = S.smiles
        SMILES = smiles
        while 1:
            hmdb = []
            colmar = []
            if SMILES in smiles_hmdb:
                hmdb = smiles_hmdb[SMILES]
                for h in hmdb:
                    if h in hmdb_to_colmar:
                        colmar.extend(hmdb_to_colmar[h])
            if SMILES in smiles_bmrb:
                colmar.extend(smiles_bmrb[SMILES])
            if colmar or SMILES == can_smiles:
                break
            SMILES = can_smiles

        hmdb = '['+', '.join(['"'+h+'"' for h in hmdb]) +']'
        colmar = '['+', '.join(['"'+h+'"' for h in colmar]) +']'

        
        motifs, spinsystems = {}, {ss: [] for ss in S.spinsystems()}
        for i in range(0, 3):
            motifs[i] = []
            motif_dict = S.motifs(i)
            for spinsys in motif_dict:
                M = Compound(motif_dict[spinsys])
                mapping = M.canonicalize()
                ss = ' '.join(sorted([str(mapping[int(n)]) for n in spinsys.split()], key=lambda x: int(x)))
                motif = M.smiles+'_'+ss
                motifs[i].extend([motif])
                spinsystems[spinsys].extend([motif])

        submotifs, nodes = {i: [] for i in range(1, 5)}, {}
        for spinsys in spinsystems:
            for node in spinsys.split():
                nodes[node] = []
                submotif_dict = S.sub_motifs(int(node))
                for shell in sorted(list(submotif_dict.keys()), key=lambda x: x):
                    SM = Compound(submotif_dict[shell])
                    mapping = SM.canonicalize()
                    submotif = SM.smiles+'_'+str(mapping[int(node)])
                    nodes[node].extend([submotif])
                    submotifs[shell].extend([submotif])
        
        spinsystems = '['+', '.join(['"'+', '.join([ss]+spinsystems[ss])+'"' for ss in spinsystems])+']'
        nodes = '['+', '.join(['"'+', '.join([n]+nodes[n])+'"' for n in nodes])+']'

        motifs = {'motif_'+str(i): '['+', '.join(['"'+m+'"' for m in motifs[i]])+']' for i in motifs}
        submotifs = {'submotif_'+str(i): '['+', '.join(['"'+m+'"' for m in submotifs[i]])+']' for i in submotifs}
        
        command = 'MATCH (m:ReferenceMolecule) WHERE m.displayName = "'+compound+'" SET m.CHEBI = "'+chebi+'", m.SMILES_3D = "'+smiles.replace('\\', '').replace('/', '')+'", m.SMILES_2D = "'+can_smiles+'", m.HMDB = '+hmdb+', m.COLMAR = '+colmar+', m.spinsystems = '+spinsystems+', m.nodes = '+nodes
        command = command+', '+', '.join(['m.'+m+' = '+motifs[m] for m in motifs])+', '+', '.join(['m.'+m+' = '+submotifs[m] for m in submotifs])
        db.run(command)


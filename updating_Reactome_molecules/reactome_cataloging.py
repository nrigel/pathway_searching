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

with open('chebi_smiles.pkl', 'rb') as pkl:
    chebi_smiles = pickle.load(pkl)


with Reactome._driver.session() as db:
    for compound in reactome_to_chebi:
        print(compound)
        chebi = reactome_to_chebi[compound]
        if chebi not in chebi_smiles:
            continue
        smiles = chebi_smiles[chebi]
        S = Compound(smiles)
        
        if not S.spinsystems(): # skip molecules that don't have carbons because they are irrelavant to us
            continue
        
        S.canonicalize()
        can_smiles = S.smiles
        
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
        
        command = 'MATCH (m:ReferenceMolecule) WHERE m.displayName = "'+compound+'" SET m.SMILES_2D = "'+smiles+'", m.SMILES_3D = "'+can_smiles+'", m.spinsystems = '+spinsystems+', m.nodes = '+nodes
        command = command+', '+', '.join(['m.'+m+' = '+motifs[m] for m in motifs])+', '+', '.join(['m.'+m+' = '+submotifs[m] for m in submotifs])
        db.run(command)
        

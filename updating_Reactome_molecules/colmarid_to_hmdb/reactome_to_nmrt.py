
dbsettings = {'NMRT': {'port': '7687', 'user': 'neo4j', 'password': 'olivia05'},
      'Reactome': {'port': '7688', 'user': 'neo4j', 'password': 'olivia05'}}
from neo4j import GraphDatabase
class NMRTServer(object):
    def __init__(self, uri='bolt://127.0.0.1:'+dbsettings['NMRT']['port']+'/', user=dbsettings['NMRT']['user'], password=dbsettings['NMRT']['password']):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))   
    def close(self):
        self._driver.close()
class ReactomeServer(object):
    def __init__(self, uri='bolt://127.0.0.1:'+dbsettings['Reactome']['port']+'/', user=dbsettings['Reactome']['user'], password=dbsettings['Reactome']['password']):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
    def close(self):
        self._driver.close()
NMRT, Reactome = NMRTServer(), ReactomeServer()

import pybel, openbabel

def remove_proton_charges(SMILES):
    # remove proton-related charges for better NMR-centric searching
    chebi = pybel.readstring("smi", SMILES)
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
    return chebi.write('can').split()[0]

with Reactome._driver.session() as reactome, NMRT._driver.session() as nmrt:
    nmrt_smiles = {}
    for node in nmrt.run('MATCH (m:metabolites) RETURN m').graph().nodes:
        smiles = None
        if 'dbIsomericSMILES' in node and node['dbIsomericSMILES'] != '':
            smiles = node['dbIsomericSMILES']
        else:
            smiles = node['dbCanonicalSMILES']
        
        nmrt_smiles[node['Compound']] = smiles

        if smiles:
            smiles = remove_proton_charges(smiles)
            if smiles not in nmrt_smiles:
                nmrt_smiles[smiles] = []
            nmrt_smiles[smiles].extend([node['Compound']])
    for node in reactome.run('MATCH (m:ReferenceMolecule) RETURN m').graph().nodes:
        #if node['COLMAR'] and len(node['COLMAR']) != len(set(node['COLMAR'])):
            #reactome.run('MATCH (m:ReferenceMolecule) WHERE m.dbId = '+str(node['dbId'])+' SET m.COLMAR = ['+', '.join(['"'+c+'"' for c in set(node['COLMAR'])])+']')

        smiles = node['SMILES_3D']
        if smiles:
            smiles = remove_proton_charges(smiles)
            if smiles in nmrt_smiles:
                colmar = []
                for compound in nmrt_smiles[smiles]:
                    if compound not in node['COLMAR']:
                        reactome.run('MATCH (m:ReferenceMolecule) WHERE m.dbId = '+str(node['dbId'])+' SET m.COLMAR = ['+', '.join(['"'+c+'"' for c in set(nmrt_smiles[smiles])])+']')
                        break

            

NMRT.close()
Reactome.close()
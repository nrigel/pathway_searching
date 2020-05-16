import xml.etree.ElementTree as ET
import networkx as nx
import csv, json, os

class PathBankDB_Builder(object):
    
    def __init__(self):
        self.node_csvs, self.edge_csvs = [], []

    def NodeWriter(self, pathbank_datadict, path):
        node_data, pathway_species  = {}, {}
        
        for node in ['Metabolite', 'Protein', 'Pathway']:
            print('Modifying '+node+' nodes...')
            node_data[node] = {}
            self.node_csvs.extend([self.path+'neo4jnodes_'+node+'.csv'])
            desc = {}

            with open(self.path+pathbank_nodes[node]['CSV']) as csvin, open(self.node_csvs[-1], 'w') as csvout:
                writer = csv.writer(csvout)
                writer.writerow(pathbank_nodes[node]['header'])
                csvin.readline()
                for row in csv.reader(csvin):
                    ID = row[pathbank_nodes[node]['id_index']]
                    write_row = row+[node]

                    if node != 'Pathway':
                        pathway_species[row[0]] = row[pathbank_nodes[node]['species_index']]
                    else:
                        if ID not in pathway_species:
                            species = 'Homo sapiens'
                        else:
                            write_row.insert(pathbank_nodes[node]['species_index'], pathway_species[ID])

                    write_row = [json.dumps(r.replace('"', '')) for r in write_row]
                    write_row = ','.join(write_row)+'\n'
                    if node == 'Pathway':
                        csvout.write(write_row)
                    else:
                        if node == 'Protein':
                            d = tuple(row[4:11])
                            if ID in node_data[node]:
                                m = len([i for i in d if i == ''])
                                M = len([i for i in desc[ID] if i == ''])
                                if m < M:
                                    desc[ID] = d
                                    node_data[node][ID] = write_row
                            else:
                                desc[ID] = d
                                node_data[node][ID] = write_row
                        else:
                            node_data[node][ID] = write_row
                if node != 'Pathway':
                    for ID in node_data[node]:
                        csvout.write(node_data[node][ID])

    def BiopaxLister(self, path):
        if path[-1] != '/':
            path = path+'/'
        filelist = []
        for PATH, dirs, files in os.walk(path):
            for filename in files:
                if filename[-3:] == 'owl':
                    filelist.extend([path+filename])
        return filelist

    def BiopaxParser(self, filename):
        level = -1
        G = nx.DiGraph() # Add data to Graph object
        skip_elements = ['Ontology', 'imports']

        for (event, elem) in ET.iterparse(filename, ['start', 'end']): 

            if event == 'end':
                level -= 1
                continue
            if event == 'start':
                level += 1
            
            if level < 1:
                continue

            TAG = elem.tag
            if '{' == TAG[0]:
                TAG = TAG.split('}')[1]
            
            if TAG in skip_elements:
                continue

            elem.set('pyindex', len(G.nodes))
            elem.set('tag', TAG)
            elem.set('level', level)
            
            if level == 1:
                parent = elem.get('pyindex')
            else:
                elem.set('XML_parent1', parent)    
            
            elem_data = {}

            if elem.text:
                elem.set('text', elem.text)

            for key, value in elem.items():
                KEY, VALUE = key, value
                if '{' == key[0]:
                    KEY = key.split('}')[1]
                
                if KEY == 'resource':
                    elem_data['type'] = value[1:].split('/')[0]
                    VALUE = VALUE[1:]

                if KEY == 'ID':
                    elem_data['type'] = value.split('/')[0]
                
                elem_data[KEY] =  VALUE

            G.add_nodes_from([(elem.get('pyindex'), elem_data)])

            if parent != elem.get('pyindex'):
                G.add_edge(parent, elem.get('pyindex'), type='XML_parent1')

        IDs = {G.nodes[node]['ID']: node for node in G.nodes if 'ID' in G.nodes[node]}

        for node in G:
            if 'resource' in G.nodes[node]:
                ref = G.nodes[node]['resource']
                if ref in IDs:
                    G.add_edge(node, IDs[ref], type='referenced')

        g = nx.DiGraph()

        pathways = [node for node in G.nodes if G.nodes[node]['tag'] == 'Pathway']

        for pathway in pathways:
            pathwayOrder = [e[1] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'pathwayOrder']
            pathwayComponent = [e[1] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'pathwayComponent']
            displayName = [G.nodes[e[1]]['text'] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'displayName']
            name = [G.nodes[e[1]]['text'] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'name']
            organism = [G.nodes[e[1]]['resource'] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'organism']
            xref = [G.nodes[e[1]]['resource'][len('Reference/SMPDB_'):] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'xref' and 'Reference/SMPDB_' in G.nodes[e[1]]['resource']]
            if not xref:
                continue
            for node in xref:
                pass

        # node IDs: Pathway -> SMPDB_ID, Protein -> Uniprot_ID, Metabolite -> PW_ID
        # pathwayOrder --> BiochemicalPathwayStep --> stepConversion --> BioChemical Reaction: left, right, name, spontaneous?, stoichiometry, conversionDirection
        # make reaction SVG images using OpenBabel, RDKit?


    def Neo4j_Command(self, neo4j_admin_path):
        neo4j_admin_command = [neo4j_admin_path+'neo4j-admin import ']
        for file in self.node_csvs:
            neo4j_admin_command.extend(['--nodes='+file])
        for file in self.edge_csvs:
            neo4j_admin_command.extend(['--relationships='+file])
        print(' '.join(neo4j_admin_command))


DBbuilder = PathBankDB_Builder()

# Node writer
pathbank_nodes = {'Pathway': {'CSV': 'pathbank_pathways.csv', 'id_index': 0, 'species_index': 4, 'header': ['SMPDB_ID:ID', 'PW_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Description', ':LABEL']},
                'Metabolite': {'CSV': 'pathbank_all_metabolites.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'PW_ID:ID', 'Metabolite_Name', 'HMDB_ID', 'KEGG_ID', 'ChEBI_ID', 'DrugBank_ID', 'CAS', 'Formula', 'IUPAC', 'SMILES', 'InChI', 'InChIKey', ':LABEL']},
                'Protein': {'CSV': 'pathbank_all_proteins.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Uniprot_ID:ID', 'Protein_Name', 'HMDBP_ID', 'DrugBank_ID', 'GenBank_ID', 'Gene_Name', 'Locus', ':LABEL']}}
#DBbuilder.NodeWriter(pathbank_datadict, path='pathbank_data/')

# List Biopax files
#biopax_list = DBbuilder.BiopaxLister(path='pathbank_data/pathbank_all_biopax/')

# Parse Biopax files
for file in ['PW002044.owl']: #biopax_list
    DBbuilder.BiopaxParser(filename=file)

import xml.etree.ElementTree as ET
import networkx as nx
import csv, json, os, traceback, types
from neo4j import GraphDatabase
import openbabel, pybel, chemstructure

class PathBankDB_Builder(object):
    
    def __init__(self):
        self.node_csvs, self.edge_csvs = [], []
        self.gene_sequences, self.protein_sequences = {}, {}
        self.pathways, self.metabolites, self.proteins = set(), set(), set()

        self.biopax_nodes, self.biopax_nodes_header = [], ['Biopax_ID:ID', 'SMPDB_ID', 'Pathway_ID', 'displayName', ':LABEL']
        self.biopax_edges, self.biopax_edges_headers = {'order': [], 'stoichiometry': [], 'else': []}, {'order': [':START_ID','order','pathway',':END_ID',':TYPE'], 'stoichiometry': [':START_ID','stoichiometry','pathway',':END_ID',':TYPE'], 'else': [':START_ID','pathway',':END_ID',':TYPE']} # ':START_ID','spinSystem',':END_ID',':TYPE'


    def SequenceReader(self, gene_fasta, protein_fasta):
        
        def fasta_reader(fasta):
            data = {}
            with open(fasta) as file:
                uniprot = None
                for line in file:
                    if '>' == line[0]:
                        uniprot = line.rstrip().split()[-1].split('(')[-1][:-1]
                        data[uniprot] = ''
                        continue
                    data[uniprot] = data[uniprot]+line.rstrip()
            return data

        for key, value in fasta_reader(gene_fasta).items():
            self.gene_sequences[key] = value
        for key, value in fasta_reader(protein_fasta).items():
            self.protein_sequences[key] = value

    def MetaboliteDescReader(self, CSV):
        result = {}
        with open(CSV, 'r') as csvin:
            for row in csv.reader(csvin):
                header = row
                break

            for row in csv.reader(csvin):
                result[row[0]] = {}
                for h in range(1, len(header)):
                    if row[h]:
                        result[row[0]][header[h]] = json.loads(row[h])
                    else:
                        result[row[0]][header[h]] = ''
        
        return result

    def NodeWriter(self, pathbank_datadict, path):
        node_data, pathway_species  = {}, {}
        nodeset = {'Metabolite': self.metabolites, 'Pathway': self.pathways, 'Protein': self.proteins}
        for node in ['Metabolite', 'Protein', 'Pathway']:
            print('Modifying '+node+' nodes...')
            node_data[node] = {}
            self.node_csvs.extend([path+'neo4j_nodes_'+node+'.csv'])
            desc = {}

            with open(path+pathbank_datadict[node]['CSV']) as csvin, open(self.node_csvs[-1], 'w') as csvout:
                writer = csv.writer(csvout)
                writer.writerow(pathbank_datadict[node]['header'])
                csvin.readline()
                for row in csv.reader(csvin):
                    ID = row[pathbank_datadict[node]['id_index']]
                    nodeset[node].add(ID)
                    write_row = row

                    if node != 'Pathway':
                        pathway_species[row[0]] = row[pathbank_datadict[node]['species_index']]
                    else:
                        if ID not in pathway_species:
                            species = 'Homo sapiens'
                        else:
                            write_row.insert(pathbank_datadict[node]['species_index'], pathway_species[ID])
                    
                    if node == 'Protein':
                        gene_seq, protein_seq = '', ''
                        if ID in self.gene_sequences:
                            gene_seq = self.gene_sequences[ID]
                        if ID in self.protein_sequences:
                            protein_seq = self.protein_sequences[ID]
                        write_row.insert(-2, gene_seq)
                        write_row.insert(-2, protein_seq)
                    
                    
                    write_row = [json.dumps(r.replace('"', '')) for r in write_row]

                    if node == 'Metabolite':
                        motiflist = {}
                        if ID in pathbank_datadict[node]['metabolite_motifs']:
                            motiflist = pathbank_datadict[node]['metabolite_motifs'][ID]
                        matchlist = {}
                        if ID in pathbank_datadict[node]['metabolite_matches']:
                            matchlist = pathbank_datadict[node]['metabolite_matches'][ID]
                        # 'spinsystems', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4', 'COLMARm'
                        for opt in ['spinsystems', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4']:
                            if opt in motiflist:
                                write_row.extend([json.dumps(';'.join(motiflist[opt]))])
                            else:
                                write_row.extend([json.dumps('')])
                        if 'COLMARm' in matchlist:
                            write_row.extend([json.dumps(';'.join(matchlist['COLMARm']))])
                            write_row.extend([json.dumps(node+';COLMARm')])
                        else:
                            write_row.extend([json.dumps('')])
                            write_row.extend([json.dumps(node)])

                    if node != 'Metabolite':
                        write_row.extend([json.dumps(node)])
                    
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

            elem.set('text', elem.text)

            for key, value in elem.items():
                KEY, VALUE = key, value
                if '{' == key[0]:
                    KEY = key.split('}')[1]
                
                if KEY == 'resource':
                    elem_data['type'] = value[1:].split('/')[0]
                    VALUE = VALUE.replace('#', '')

                if KEY == 'ID':
                    elem_data['type'] = value.split('/')[0]
                
                elem_data[KEY] =  VALUE

            if TAG == 'Pathway':
                if 'about' in elem_data and 'ID' not in elem_data:
                    elem_data['ID'] = elem_data['about'] # pathway references are based on the 'about' element

            G.add_nodes_from([(elem.get('pyindex'), elem_data)])

            if parent != elem.get('pyindex'):
                G.add_edge(parent, elem.get('pyindex'), type='XML_parent1')

        IDs = {G.nodes[node]['ID']: node for node in G.nodes if 'ID' in G.nodes[node]}
        # merge resources with sources and then remove resource nodes
        for node in list(G.nodes):
            if 'resource' in G.nodes[node]:
                ref = G.nodes[node]['resource']
                if ref in IDs:
                    pair = (list(G.in_edges(node))[0][0], IDs[ref])
                    node_dict = {pair: dict(G.nodes[node])}
                    # 'IDs[ref]' is the actual node, 'node' is the reference
                    # take all in_edge to 'node' and make them with 'IDs[ref]'
                    # make out_edge
                    G.add_edge(pair[0], pair[1])
                    nx.set_edge_attributes(G, node_dict)
                    G.remove_node(node)

        g = nx.DiGraph()
        file_id = filename.split('/')[-1].split('.')[0]
        master_pathway = None
        # compile graph
        for pathway in [node for node in G.nodes if G.nodes[node]['tag'] == 'Pathway']:
            n = self.Biopax_Component(Pathway=pathway, G=G)
            g.add_nodes_from([n])
            pathwayID = n[0]

            xrefs = {G.nodes[e[1]]['ID'][len('Reference/'):].split('_')[0]: G.nodes[e[1]]['ID'][len('Reference/'):].split('_')[1] for e in G.out_edges(pathway) if 'tag' in G.nodes[e[1]] and G.nodes[e[1]]['tag'] == 'UnificationXref'}
            if 'PathWhiz' in xrefs:
                if xrefs['PathWhiz'] == file_id:
                    master_pathway = pathwayID
            
            pathwayOrder, pathwayComponent = sorted([(e[1], G.edges[e]['pyindex']) for e in G.out_edges(pathway) if 'tag' in G.edges[e] and G.edges[e]['tag'] == 'pathwayOrder'], key=lambda x: x[1]), sorted([(e[1], G.edges[e]['pyindex']) for e in G.out_edges(pathway) if 'tag' in G.edges[e] and G.edges[e]['tag'] == 'pathwayComponent'], key=lambda x: x[1])
            prev_node, order = None, 0
            for node, ref in pathwayOrder:
                for e1 in G.out_edges(node):
                    e1_nodedict = {'type': G.nodes[e1[1]]['tag'], 'displayName': ''}

                    # Add in metabolites from the reaction

                    if G.nodes[e1[1]]['tag'] in ['BiochemicalReaction', 'Transport', 'TransportWithBiochemicalReaction']:

                        for e2 in G.out_edges(e1[1]): 
                            if G.nodes[e2[1]]['tag'] in ['spontaneous', 'conversionDirection']:
                                e1_nodedict[G.nodes[e2[1]]['tag']] = G.nodes[e2[1]]['text']
                        
                        if G.nodes[e1[1]]['ID'] not in g.nodes:
                            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(e1[1]) if G.nodes[e2[1]]['tag'] == 'name']
                            if name:
                                e1_nodedict['displayName'] = name[0]
                            g.add_nodes_from([(G.nodes[e1[1]]['ID'], e1_nodedict)])
                        
                        order += 1
                        g.add_edge(pathwayID, G.nodes[e1[1]]['ID'], type='hasEvent', order=order)

                        if prev_node:
                           g.add_edge(G.nodes[e1[1]]['ID'], prev_node, type='preceedingEvent') 
                        prev_node = G.nodes[e1[1]]['ID']

                        conversion = {'left': 'input', 'right': 'output'}
                        if 'conversionDirection' in e1_nodedict and e1_nodedict['conversionDirection'] == 'RIGHT_TO_LEFT':
                            conversion = {'left': 'output', 'right': 'input'}

                        stoichiometry = {}
                        for stoich in [e2[1] for e2 in G.out_edges(e1[1]) if G.nodes[e2[1]]['tag'] == 'Stoichiometry']:
                            stoichiometricCoefficient = [G.nodes[e3[1]]['text'] for e3 in G.out_edges(stoich) if G.nodes[e3[1]]['tag'] == 'stoichiometricCoefficient'][0]
                            molecule = [G.nodes[e3[1]]['ID'] for e3 in G.out_edges(stoich) if G.nodes[e3[1]]['tag'] != 'stoichiometricCoefficient'][0]
                            stoichiometry[molecule] = stoichiometricCoefficient

                        for e2 in G.out_edges(e1[1]):
                            if G.nodes[e2[1]]['tag'] == 'SmallMolecule':
                                n = self.Biopax_Component(Compound=e2[1], G=G)
                                if not n:
                                    continue
                                if n[0] not in g.nodes:
                                    g.add_nodes_from([n])
                                
                                edge_dict = {'type': conversion[G.edges[e2]['tag']], 'stoichiometry': stoichiometry[G.nodes[e2[1]]['ID']]}
                                
                                if edge_dict['type'] == 'input': # e2[1] -> e1[1]
                                    pair = (n[0], G.nodes[e1[1]]['ID'])
                                if edge_dict['type'] == 'output': # e1[1] -> e2[1]
                                    pair = (G.nodes[e1[1]]['ID'], n[0])
                                
                                g.add_edge(pair[0], pair[1])
                                nx.set_edge_attributes(g, {pair: edge_dict})
                        
                    # add in proteins from the reaction
                    elif G.nodes[e1[1]]['tag'] == 'Catalysis':

                        e1_nodedict['controlled'] = [G.nodes[e2[1]]['ID'] for e2 in G.out_edges(e1[1]) if 'tag' in G.edges[e2] and G.edges[e2]['tag'] == 'controlled'][0]
                        e1_nodedict['controlType'] = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(e1[1]) if G.nodes[e2[1]]['tag'] == 'controlType'][0]

                        if G.nodes[e1[1]]['ID'] not in g.nodes:
                            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(e1[1]) if G.nodes[e2[1]]['tag'] == 'displayName']
                            if name:
                                e1_nodedict['displayName'] = name[0]
                            g.add_nodes_from([(G.nodes[e1[1]]['ID'], e1_nodedict)])
                        
                        g.add_edge(pathwayID, G.nodes[e1[1]]['ID'], type='hasEvent', order=order)

                        if prev_node:
                           g.add_edge(G.nodes[e1[1]]['ID'], prev_node, type='catalysisEvent') 

                        controller = [e2[1] for e2 in G.out_edges(e1[1]) if 'tag' in G.edges[e2] and G.edges[e2]['tag'] == 'controller'][0]

                        if 'Protein/' in G.nodes[controller]['ID']:
                            n = self.Biopax_Component(Protein=controller, G=G)
                            if n:
                                if n[0] not in g.nodes:
                                    g.add_nodes_from([n])
                                g.add_edge(G.nodes[e1[1]]['ID'], n[0], type='hasEvent')
                            
                        elif 'ProteinComplex/' in G.nodes[controller]['ID']:
                            nodes, edges = self.Biopax_Complex(Complex=controller, G=G)
                            for n in nodes:
                                if n[0] not in g:
                                    g.add_nodes_from([n])
                            for e in edges:
                                E = list(e)
                                E[E.index('parent')] = G.nodes[e1[1]]['ID']
                                g.add_edge(E[0], E[1])
                                nx.set_edge_attributes(g, {tuple(E): edges[e]})
                        else:
                            raise Exception(G.nodes[controller]['ID'], 'catalysis issue')

                    elif G.nodes[e1[1]]['tag'] == 'Interaction' or G.nodes[e1[1]]['tag'] == 'MolecularInteraction':
                        
                        if G.nodes[e1[1]]['ID'] not in g.nodes:
                            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(e1[1]) if G.nodes[e2[1]]['tag'] == 'name']
                            if name:
                                e1_nodedict['displayName'] = name[0]
                            g.add_nodes_from([(G.nodes[e1[1]]['ID'], e1_nodedict)])
                        
                        order += 1
                        g.add_edge(pathwayID, G.nodes[e1[1]]['ID'], type='hasEvent', order=order)

                        if prev_node:
                           g.add_edge(G.nodes[e1[1]]['ID'], prev_node, type='preceedingEvent') 
                        prev_node = G.nodes[e1[1]]['ID']

                        for participant in [e2[1] for e2 in G.out_edges(e1[1]) if 'tag' in G.edges[e2] and G.edges[e2]['tag'] == 'participant']:
                            if G.nodes[participant]['tag'] == 'SmallMolecule':
                                n = self.Biopax_Component(Compound=participant, G=G)
                                if not n:
                                    continue
                                if n[0] not in g.nodes:
                                    g.add_nodes_from([n])
                                g.add_edge(G.nodes[e1[1]]['ID'], n[0], type='hasEvent')
                            
                            elif G.nodes[participant]['tag'] == 'Complex':
                                nodes, edges = self.Biopax_Complex(Complex=participant, G=G)
                                for n in nodes:
                                    if n[0] not in g:
                                        g.add_nodes_from([n])
                                for e in edges:
                                    E = list(e)
                                    E[E.index('parent')] = G.nodes[e1[1]]['ID']
                                    g.add_edge(E[0], E[1])
                                    nx.set_edge_attributes(g, {tuple(E): edges[e]})

                            elif G.nodes[participant]['tag'] == 'Protein':
                                n = self.Biopax_Component(Protein=participant, G=G)
                                if n:
                                    if n[0] not in g.nodes:
                                        g.add_nodes_from([n])
                                    g.add_edge(G.nodes[e1[1]]['ID'], n[0], type='hasEvent')

                            elif G.nodes[participant]['tag'] == 'Pathway':
                                n = self.Biopax_Component(Pathway=participant, G=G)
                                if n:
                                    if n[0] not in g.nodes:
                                        g.add_nodes_from([n])
                                    g.add_edge(G.nodes[e1[1]]['ID'], n[0], type='hasEvent')
                            
                            elif G.nodes[participant]['tag'] in ['Dna', 'Rna']:
                                n = self.Biopax_Component(NA=participant, G=G)
                                if n:
                                    if n[0] not in g.nodes:
                                        g.add_nodes_from([n])
                                    g.add_edge(G.nodes[e1[1]]['ID'], n[0], type='hasEvent')
                            
                            else:
                                raise Exception(G.nodes[participant]['tag'], G.nodes[e1[1]])

                    elif G.nodes[e1[1]]['tag'] == 'Pathway':
                        n = self.Biopax_Component(Pathway=e1[1], G=G)
                        if n[0] not in g:
                            g.add_nodes_from([n])
                        
                        order += 1
                        g.add_edge(pathwayID, n[0], type='hasEvent', order=order)
                        if prev_node:
                           g.add_edge(n[0], prev_node, type='preceedingEvent') 
                        prev_node = n[0]

                    else:
                        raise Exception(G.nodes[e1[1]]['tag'])

        # clean up graph
        # nodes already present: pathways, metabolites, proteins
        # nodes specific to the pathway to be written: reactions, catalysts, transports, interactions, etc.
        nodeset = {'Metabolite': self.metabolites, 'Pathway': self.pathways, 'Protein': self.proteins}
        if master_pathway not in nodeset['Pathway']:
            print(filename)

        pathway_nodes = [node for node in g.nodes if g.nodes[node]['type'] not in ['Pathway', 'Protein', 'Metabolite']]
        for node in [node for node in g.nodes if g.nodes[node]['type'] in ['Pathway', 'Protein', 'Metabolite']]:
            if node not in nodeset[g.nodes[node]['type']]:
                pathway_nodes.extend([node])
            else:
                if g.nodes[node]['type'] == 'Metabolite':
                    self.biopax_edges['else'].extend([[master_pathway, master_pathway, node, 'hasMetabolite']])

        pathway_nodes = {pathway_nodes[i]: file_id+'-'+master_pathway+'_'+str(i+1).zfill(len(str(len(pathway_nodes)+1))) for i in range(len(pathway_nodes))}
        for node in pathway_nodes:
            displayName = ''
            if 'displayName' in g.nodes[node]:
                displayName = g.nodes[node]['displayName']
            self.biopax_nodes.extend([[ pathway_nodes[node], master_pathway, node.replace('/', '-'), displayName, g.nodes[node]['type'] ]])

        for edge in g.edges:
            node1, node2 = edge
            if node1 in pathway_nodes:
                node1 = pathway_nodes[node1]
            if node2 in pathway_nodes:
                node2 = pathway_nodes[node2]
            key = [k for k in dict(g.edges[edge]).keys() if k != 'type']
            if key:
                self.biopax_edges[key[0]].extend([[node1, g.edges[edge][key[0]], master_pathway, node2, g.edges[edge]['type']]])
            else:
                self.biopax_edges['else'].extend([[node1, master_pathway, node2, g.edges[edge]['type']]])
        
        return g          

    def Biopax_Component(self, G, Pathway=None, Compound=None, Protein=None, NA=None):
        if Pathway is not None:
            ID = G.nodes[Pathway]['ID'].split('/')[-1]
            if ID[:3] =='SMP':
                ID = 'SMP'+ID[3:].zfill(7)
            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(Pathway) if G.nodes[e2[1]]['tag'] == 'displayName']
            if name:
                name = name[0]
            else:
                name = ''
            return (ID, {'type': 'Pathway', 'displayName': name})
        
        if Compound is not None:
            ID = G.nodes[Compound]['ID'][len('Compound/'):] # get DB ID
            if 'Generic' in ID:
                return # don't think we need to include these
            if len(ID.split('_')) > 2:
                ID = '_'.join(ID.split('_')[:2])
            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(Compound) if G.nodes[e2[1]]['tag'] == 'displayName']
            if name:
                name = name[0]
            else:
                name = ''
            return (ID, {'type': 'Metabolite', 'displayName': name})
        
        if Protein is not None:
            ProteinReference = [e2[1] for e2 in G.out_edges(Protein) if 'type' in G.nodes[e2[1]] and G.nodes[e2[1]]['type'] == 'ProteinReference']
            if not ProteinReference:
                return
            ProteinReference = ProteinReference[0]
            ID = [G.nodes[e2[1]]['ID'][len('Reference/UniProt_'):] for e2 in G.out_edges(ProteinReference) if G.nodes[e2[1]]['tag'] == 'UnificationXref']
            if not ID:
                return
            ID = ID[0]
            name = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(Protein) if G.nodes[e2[1]]['tag'] == 'displayName']
            if name:
                name = name[0]
            else:
                name = ''
            return (ID, {'type': 'Protein', 'displayName': name})

        if NA is not None:
            tag = G.nodes[NA]['tag'][0]+'NA'
            displayName = [G.nodes[e2[1]]['text'] for e2 in G.out_edges(NA) if G.nodes[e2[1]]['tag'] == 'displayName']
            nareference = [e2[1] for e2 in G.out_edges(NA) if G.nodes[e2[1]]['tag'] == G.nodes[NA]['tag']+'Reference']
            db = None
            for ref in nareference:
                xref = [e2[1] for e2 in G.out_edges(ref) if G.nodes[e2[1]]['tag'] == 'UnificationXref']
                if xref:
                    db = [G.nodes[e3[1]]['text'] for e3 in G.out_edges(xref[0]) if G.nodes[e3[1]]['tag'] == 'db']
                    if db:
                        db = db[0]
                        break
            if db:
                return (displayName[0], {'type': tag, 'db': db, 'displayname': displayName[0]})

            if displayName:
                return (displayName[0], {'type': tag, 'displayname': displayName[0]})


    def Biopax_Complex(self, Complex, G):
        nodes, edges = [], {}
        stoichiometry = {G.nodes[stoich]['ID'][:G.nodes[stoich]['ID'].find('/Stoichiometry/')]: G.nodes[stoich]['ID'][G.nodes[stoich]['ID'].find('/Stoichiometry/')+len('/Stoichiometry/'):] for stoich in [e2[1] for e2 in G.out_edges(Complex) if G.nodes[e2[1]]['tag'] == 'Stoichiometry']}
        for protein in [e2[1] for e2 in G.out_edges(Complex) if G.nodes[e2[1]]['tag'] == 'Protein']:
            n = self.Biopax_Component(Protein=protein, G=G)
            if not n:
                continue
            nodes.extend([n])
            node_dict = {'type': 'hasEvent', 'stoichiometry': '1.0'}
            if G.nodes[protein]['ID'] in stoichiometry:
                node_dict['stoichiometry'] = stoichiometry[G.nodes[protein]['ID']]
            edges[('parent', nodes[-1][0])] = node_dict

        for compound in [e2[1] for e2 in G.out_edges(Complex) if G.nodes[e2[1]]['tag'] == 'SmallMolecule']:
            n = self.Biopax_Component(Compound=compound, G=G)
            if not n:
                continue
            nodes.extend([n])
            node_dict = {'type': 'hasEvent', 'stoichiometry': '1.0'}
            if G.nodes[compound]['ID'] in stoichiometry:
                node_dict['stoichiometry'] = stoichiometry[G.nodes[compound]['ID']]
            edges[('parent', nodes[-1][0])] = node_dict

        return nodes, edges

        # node IDs: Pathway -> SMPDB_ID, Protein -> Uniprot_ID, Metabolite -> PW_ID
        # pathwayOrder --> BiochemicalPathwayStep --> stepConversion --> BioChemical Reaction: left, right, name, spontaneous?, stoichiometry, conversionDirection
        # make reaction SVG images using OpenBabel, RDKit?

    def BiopaxNodeWriter(self, path):
        if self.biopax_nodes_header:
            with open(path+'neo4j_biopax_nodes.csv', 'w') as nodes:
                writer = csv.writer(nodes)
                writer.writerow(self.biopax_nodes_header)
            self.node_csvs.extend([path+'neo4j_biopax_nodes.csv'])
            self.biopax_nodes_header = []
        
        with open(path+'neo4j_biopax_nodes.csv', 'a') as nodes:
            writer = csv.writer(nodes)
            for node in self.biopax_nodes:
                writer.writerow(node)
            self.biopax_nodes = []

    def BiopaxEdgeWriter(self, path):
        for edge in list(self.biopax_edges_headers.keys()):
            with open(path+'neo4j_biopax_edges_'+edge+'.csv', 'w') as edges:
                writer =csv.writer(edges)
                writer.writerow(self.biopax_edges_headers[edge])
            self.edge_csvs.extend([path+'neo4j_biopax_edges_'+edge+'.csv'])
            del self.biopax_edges_headers[edge]
        
        for edge in self.biopax_edges:
            with open(path+'neo4j_biopax_edges_'+edge+'.csv', 'a') as edges:
                writer =csv.writer(edges)
                for e in self.biopax_edges[edge]:
                    writer.writerow(e)
                self.biopax_edges[edge] = []


    def Neo4j_Command(self, neo4j_admin_path):
        if neo4j_admin_path[-1] != '/':
            neo4j_admin_path = neo4j_admin_path+'/'
        neo4j_admin_command = [neo4j_admin_path+'bin/neo4j-admin import ']
        for file in self.node_csvs:
            neo4j_admin_command.extend(['--nodes='+file])
        for file in self.edge_csvs:
            neo4j_admin_command.extend(['--relationships='+file])
        print(' '.join(neo4j_admin_command))

    def MetaboliteMatching(self, NMRT, PathBank):
        
        def remove_charge(SMILES):
            mol = pybel.readstring("smi", SMILES)
            can = mol.write('can')
            if mol.OBMol.GetTotalCharge() or '+' in can or '-' in can:
                mol.OBMol.BeginModify()
                mol.OBMol.SetTotalCharge(0)
                for atom in openbabel.OBMolAtomIter(mol.OBMol):
                    charge = atom.GetFormalCharge()
                    if charge > 0:
                        if atom.ImplicitHydrogenCount():
                            atom.SetFormalCharge(0)
                    if charge < 0:
                        atom.SetFormalCharge(0)
                mol.OBMol.EndModify()
            can = mol.write('can').split()
            if can:
                return can[0]

        def remove_stereochemistry(SMILES):
            mol = pybel.readstring("smi", SMILES)

            mol.make3D() # generate 3D coordinates to avoid error message printed in terminal
            mol.removeh() # remove hydrogens
            mol.draw(show=False, filename=None, update=True, usecoords=False) # generate 2D coordinates
            can = mol.write('can').split()
            if can:
                return can[0]


        with NMRT._driver.session() as nmrt, PathBank._driver.session() as pathbank:
            metabolites = {'NMRT': nmrt.run('MATCH (m:metabolites) RETURN ([m.Compound, m.SMILES])').value(), 
                        'PathBank': [metabolite for metabolite in pathbank.run('MATCH (m:Metabolite) RETURN ([m.PW_ID, m.SMILES])').value() if metabolite[0]]}
            
            nmrt_smiles = {}
            for metabolite in metabolites['NMRT']:
                smiles = remove_charge(metabolite[1])
                if smiles not in nmrt_smiles:
                    nmrt_smiles[smiles] = []
                nmrt_smiles[smiles].extend([metabolite[0]])
            
            for metabolite in metabolites['PathBank']:
                if metabolite[1] is None or metabolite[1] == '':
                    continue
                if 'c' not in metabolite[1] and 'C' not in metabolite[1]:
                    continue
                
                structure = chemstructure.Compound(metabolite[1])
                structure.canonicalize()
                smiles = structure.smiles
                if smiles in nmrt_smiles:
                    pathbank.run('MATCH (m:Metabolite) WHERE m.PW_ID = "'+metabolite[0]+'" SET m.COLMARm = '+json.dumps(nmrt_smiles[smiles])+', m :COLMARm')
                
    def AddMotifs(self, PathBank, path, Restart=True):
        if path[-1] != '/':
            path = path+'/'

        SEEN = []
        if Restart is False:
            with open(path+'pathbank_motiflists.csv', 'w') as csvout:
                writer = csv.writer(csvout)
                writer.writerow(['PW_ID', 'spinsystems', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4'])
        else:
            with open(path+'pathbank_motiflists.csv', 'r') as csvin:
                csvin.readline()
                for row in csv.reader(csvin):
                    SEEN.extend([row[0]])

        with PathBank._driver.session() as pathbank:
            metabolites = [metabolite for metabolite in pathbank.run('MATCH (m:Metabolite) RETURN ([m.PW_ID, m.SMILES])').value() if metabolite[0] and metabolite[0] not in SEEN]

            for metabolite in metabolites:
                with open(path+'pathbank_motiflists.csv', 'a') as csvout:
                    writer = csv.writer(csvout)
                    
                    if metabolite[1] is None or metabolite[1] == '':
                        writer.writerow([metabolite[0], '', '', '', '', '', '', '', ''])
                        continue
                    if 'c' not in metabolite[1] and 'C' not in metabolite[1]:
                        writer.writerow([metabolite[0], '', '', '', '', '', '', '', ''])
                        continue
                    
                    structure = chemstructure.Compound(metabolite[1])
                    structure.canonicalize()
                    
                    motiflists = [{}, {}] # by shell, by spinsystem
                    submotiflists = {1: set(), 2: set(), 3: set(), 4: set()} # by shell
                    for shell in [0, 1, 2]:
                        motiflists[0][shell] = set()
                        motifs = structure.motifs(shell)
                        for spinsys in motifs:
                            if spinsys not in motiflists[1]:
                                motiflists[1][spinsys] = []

                            S = chemstructure.Compound(motifs[spinsys])
                            aliases = S.canonicalize()
                            new_spinsys = ' '.join(sorted([str(aliases[int(n)]) for n in spinsys.split()], key=lambda x: int(x)))
                            motif = S.smiles+'_'+new_spinsys

                            motiflists[0][shell].add(motif)
                            motiflists[1][spinsys].extend([motif])
                    motiflists[0] = {shell: json.dumps(list(motiflists[0][shell])) for shell in motiflists[0]}
                    motiflists[1] = json.dumps(['.'.join([spinsys]+motiflists[1][spinsys]) for spinsys in motiflists[1]])

                    for spinsys in structure.spinsystems():
                        for node in spinsys.split():
                            submotifs = structure.sub_motif_smiles(int(node))
                            for shell in submotiflists:
                                submotiflists[shell].add(submotifs[shell])

                    submotiflists = {shell: json.dumps(list(submotiflists[shell])) for shell in submotiflists}

                    write_list = [metabolite[0], motiflists[1]]
                    for shell in [0, 1, 2]:
                        write_list.extend([motiflists[0][shell]])
                    for shell in [1, 2, 3, 4]:
                        write_list.extend([submotiflists[shell]])

                    writer.writerow(write_list)

    def AddSVG(self, PathBank):
        pass

    def AddReactionLabel(self, PathBank):
        with PathBank._driver.session() as db:
            for node in ['BiochemicalReaction', 'TransportWithBiochemicalReaction', 'Transport', 'Interaction', 'MolecularInteraction']:
                db.run('MATCH (n:'+node+') SET n :Reaction')


class Neo4jServer(object):
    def __init__(self, dbsettings):
        uri = 'bolt://127.0.0.1:'+dbsettings['port']+'/'
        user, password = dbsettings['user'], dbsettings['password']
        self._driver = GraphDatabase.driver(uri, auth=(user, password))   
    def close(self):
        self._driver.close()



DBbuilder = PathBankDB_Builder()

CURATE = False
if CURATE:

    # Read in sequences
    gene_fasta = '/Users/nickrigel/Documents/GitHub/pathway_searching/curate_database/pathbank_data/pathbank_gene.fasta'
    protein_fasta = '/Users/nickrigel/Documents/GitHub/pathway_searching/curate_database/pathbank_data/pathbank_protein.fasta'
    DBbuilder.SequenceReader(gene_fasta, protein_fasta)

    metabolite_motifs = DBbuilder.MetaboliteDescReader('/Users/nickrigel/Documents/GitHub/pathway_searching/curate_database/pathbank_data/pathbank_motiflists.csv')
    metabolite_matches = DBbuilder.MetaboliteDescReader('/Users/nickrigel/Documents/GitHub/pathway_searching/curate_database/pathbank_data/pathbank_colmarm_matches.csv')

    # Node writer
    pathbank_datadict = {'Pathway': {'CSV': 'pathbank_pathways.csv', 'id_index': 0, 'species_index': 4, 'header': ['SMPDB_ID:ID', 'PW_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Description', ':LABEL']},
                    'Metabolite': {'metabolite_motifs': metabolite_motifs, 'metabolite_matches': metabolite_matches, 'CSV': 'pathbank_all_metabolites.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'PW_ID:ID', 'Metabolite_Name', 'HMDB_ID', 'KEGG_ID', 'ChEBI_ID', 'DrugBank_ID', 'CAS', 'Formula', 'IUPAC', 'SMILES', 'InChI', 'InChIKey', 'spinsystems:string[]', 'motif_0:string[]', 'motif_1:string[]', 'motif_2:string[]', 'submotif_1:string[]', 'submotif_2:string[]', 'submotif_3:string[]', 'submotif_4:string[]', 'COLMARm:string[]', ':LABEL']},
                    'Protein': {'CSV': 'pathbank_all_proteins.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Uniprot_ID:ID', 'Protein_Name', 'HMDBP_ID', 'DrugBank_ID', 'GenBank_ID', 'Gene_Name', 'Locus', 'Gene_Sequence', 'Protein_Sequence', ':LABEL']}}
    DBbuilder.NodeWriter(pathbank_datadict, path='pathbank_data/')

    DBbuilder.Neo4j_Command(neo4j_admin_path='/Applications/neo4j-community-3.5.14_pathbank/')
    exit()
    # List Biopax files
    biopax_list = DBbuilder.BiopaxLister(path='pathbank_data/pathbank_all_biopax/')

    # Parse Biopax files
    for f in range(len(biopax_list)): #['PW002044.owl']
        if biopax_list[f] != 'pathbank_data/pathbank_all_biopax/PW000965.owl':
            #continue
            pass
        try:
            DBbuilder.BiopaxParser(filename=biopax_list[f])
            DBbuilder.BiopaxNodeWriter(path='pathbank_data/')
            DBbuilder.BiopaxEdgeWriter(path='pathbank_data/')
        except Exception as e:
            print(traceback.format_exc())
            print(f+1, biopax_list[f])
            exit()
    DBbuilder.Neo4j_Command(neo4j_admin_path='/Applications/neo4j-community-3.5.14_pathbank/')

    # pathbank_gene.fasta -> uniprot to RNA sequence
    # pathbank_protein.fasta -> uniprot to amino acid sequence
    # pathbank_all_metabolites.csv -> metabolite nodes
    # pathbank_all_proteins.csv -> protein nodes
    # pathbank_pathways.csv -> pathway nodes
    # all_pathway_reactions -> reaction structures for each species

PROCESS = True
if PROCESS:
    dbsettings = {'NMRT': {'port': '7687', 'user': 'neo4j', 'password': 'olivia05'},
                'PathBank': {'port': '7690', 'user': 'neo4j', 'password': 'olivia05'}}

    NMRT = Neo4jServer(dbsettings['NMRT'])
    PathBank = Neo4jServer(dbsettings['PathBank'])

    # Match metabolites to NMRT
    #DBbuilder.MetaboliteMatching(NMRT, PathBank)

    # Add motifs and sub-motifs to metabolites
    #DBbuilder.AddMotifs(PathBank, path='pathbank_data/')

    # Add label to reaction-type nodes
    DBbuilder.AddReactionLabel(PathBank)

    # 255 Reactome matches -> 436 PathBank matches


# Fix removal of slashes from SMILES in the node writer function


#/Applications/neo4j-community-3.5.14_pathbank/bin/neo4j-admin import  --nodes=pathbank_data/neo4j_nodes_Metabolite.csv --nodes=pathbank_data/neo4j_nodes_Protein.csv --nodes=pathbank_data/neo4j_nodes_Pathway.csv --nodes=pathbank_data/neo4j/neo4j_biopax_nodes.csv --relationships=pathbank_data/neo4j/neo4j_biopax_edges_else.csv --relationships=pathbank_data/neo4j/neo4j_biopax_edges_order.csv --relationships=pathbank_data/neo4j/neo4j_biopax_edges_stoichiometry.csv
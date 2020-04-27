#!/opt/anaconda2/envs/P3/bin/python
import cgi, cgitb, sys, math, csv, json
import networkx as nx
sys.path.append('../../scripts')
from neo4j_connect import NMRTServer, ReactomeServer

# Create instance of FieldStorage 
form = cgi.FieldStorage() 

print("Content-type:text/html\r\n\r\n")

# Get data from fields
datafile = form['userfile'].file
structure_opts = [o for o in ['metabolites', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4'] if form.getvalue(o)]
species = form.getvalue('specieslist')
count_cutoff = form.getvalue('cutoff_count')

# Read in data file
lines = []
while True:
    line = datafile.readline()
    if not line:
        break
    lines.extend([line])
data_set = [line.rstrip().decode('utf-8').split(',') for line in lines]

# Format metabolites dict
names_to_colmarid = {}
with open('/Users/nickrigel/Documents/GitHub/pathway_searching/updating_Reactome_molecules/name_to_colmarid.csv') as csvin:
    csvin.readline()
    for row in csv.reader(csvin):
        names_to_colmarid[row[0]] = row[1]

metabolites = {}
for i in range(1, len(data_set)):
    metabolite, p_value, fold_change = data_set[i][0], -math.log(float(data_set[i][1])), float(data_set[i][2])
    if metabolite in names_to_colmarid:
        colmar = names_to_colmarid[metabolite]
    else:
        colmar = metabolite    
    metabolites[colmar] = {'COLMAR': colmar, 'Fold Change': fold_change, 'P-Value': p_value, 'Metabolite': metabolite}

# 1. Search for pathways with any of the inputted metabolites
Reactome = ReactomeServer()
pathways = [{}, [], [], []]
with Reactome._driver.session() as db:
    for metabolite in metabolites:
        pathways[0][metabolite] = set()
        colmar = metabolites[metabolite]['COLMAR']
        fold_change = metabolites[metabolite]['Fold Change']
        for (m, p) in db.run('MATCH (m:ReferenceMolecule)-[me:referenceEntity]-(e:PhysicalEntity)-[er:input]-(r:Reaction)-[rp:hasEvent]-(p:Pathway) WHERE "'+colmar+'" IN m.COLMAR AND p.speciesName = "Homo sapiens" RETURN ([m.displayName, p.displayName])').value():
            pathways[0][metabolite].add(p)

# 2. Generate list of pathways
for metabolite in pathways[0]:
    for pathway in pathways[0][metabolite]:
        pathways[1].extend([pathway])

# 3. Sorte pathway list by most occuring coming first
pathways[1] = sorted(pathways[1], key=lambda x: -pathways[1].count(x))

# 4. Generate list of pathways in the order from step 3 with no repeats
for p in pathways[1]:
    if p not in pathways[2]:
        pathways[2].extend([p])

# 5. Iterate thru pathway list and compile new list containing pathways that have "n" number of metabolites
for p in pathways[2]:
    count = pathways[1].count(p)
    if count >= int(count_cutoff):
        pathways[3].extend([p])

# 6. Search for all reactions in each of those pathways
pathway_graphs = {}
with Reactome._driver.session() as db:
    # https://neo4j.com/docs/api/python-driver/current/types/graph.html
    entities_to_metabolites = {}
    #pathways[3] = ["Citric acid cycle (TCA cycle)"]
    pathway_names, pathway_data = [], {}
    for pathway in pathways[3]:
        pathway_graphs[pathway] = nx.DiGraph()

        command = 'MATCH (p:Pathway) WHERE p.displayName = "'+pathway+'" AND p.speciesName = "'+species+'"'
        command = command+'\n'+'MATCH (R1:Reaction)-[pe:precedingEvent]-(R2:Reaction) WHERE (R1)-[:hasEvent]-(p)  AND (R2)-[:hasEvent]-(p) '
        command = command+'\n'+'RETURN R1, R2, pe'

        reactions = db.run(command).graph()

        if reactions.nodes:
            #print(pathway, '\n')
            for node in reactions.nodes:
                if 'Reaction' in list(node.labels): 
                    node_dict = dict(node)
                    node_dict['labels'], node_dict['TYPE'] = list(node.labels), 'Reaction'
                    pathway_graphs[pathway].add_nodes_from([(node['dbId'], node_dict)])

                    command = 'MATCH (p:Pathway) WHERE p.displayName = "'+pathway+'" AND p.speciesName = "'+species+'"'
                    command = command+'\n'+'MATCH (m:ReferenceMolecule)-[me:referenceEntity]-(e:PhysicalEntity)-[er]-(r:Reaction)-[rp:hasEvent]-(p) WHERE r.dbId = '+str(node['dbId'])+' RETURN *'
                    reaction = db.run(command).graph()
                    
                    # Reaction --> PhysicalEntity --> ReferenceMolecule
                    # make dictionary of PhysicalEntities to Metabolites
                    for relationship in reaction.relationships:
                        if relationship.type == 'referenceEntity':
                            n1, n2 = relationship.start_node, relationship.end_node # n1 = entity, n2 = metabolite
                            entities_to_metabolites[n1['dbId']] = n2 

                    # add metabolites to pathway graph by connecting to reactions
                    for relationship in reaction.relationships:
                        n1, n2 = relationship.start_node, relationship.end_node # n1 = Reaction, n2 = PhysicalEntity
                        if relationship.type in ['input', 'output'] and 'Reaction' in list(n1.labels):
                            n3 = entities_to_metabolites[n2['dbId']] # metabolite node
                            if n3['dbId'] not in pathway_graphs[pathway]: # add metabolite node to new graph
                                node_dict = dict(n3)
                                node_dict['labels'], node_dict['TYPE'] = list(n3.labels), 'Metabolite'
                                pathway_graphs[pathway].add_nodes_from([(n3['dbId'], node_dict)])
                            if relationship.type == 'input': # input metabolites go towards reaction
                                pathway_graphs[pathway].add_edge(n3['dbId'], n1['dbId'], relationship=relationship, TYPE='input')
                            if relationship.type == 'output': # output metabolites go away from reaction
                                pathway_graphs[pathway].add_edge(n1['dbId'], n3['dbId'], relationship=relationship, TYPE='output')

            for relationship in reactions.relationships:
                n2, n1 = relationship.start_node, relationship.end_node
                # (n1)-[:PrecedingEvent]->(n2); start_node = reaction 2, end_node = reaction 1
                pathway_graphs[pathway].add_edge(n1['dbId'], n2['dbId'], TYPE='precedingEvent', relationship=relationship)

        # input metabolites go towards reaction
        # output metabolites go away from reaction
        # TYPE is Reaction or Metabolite

        G = nx.DiGraph()
        for node in pathway_graphs[pathway].nodes:
            if pathway_graphs[pathway].nodes[node]['TYPE'] == 'Metabolite':
                node_dict = pathway_graphs[pathway].nodes[node]
                node_dict['id'] = node
                for colmar in node_dict['COLMAR']:
                    if colmar in metabolites:
                        node_dict['Fold Change'], node_dict['P-Value'] = metabolites[colmar]['Fold Change'], metabolites[colmar]['P-Value']
                G.add_nodes_from([(node, node_dict)])
                #G.add_node(node, id=node, name=pathway_graphs[pathway].nodes[node]['displayName'], type=pathway_graphs[pathway].nodes[node]['TYPE'])
        
        for edge in pathway_graphs[pathway].edges:
            reaction = [n for n in range(len(edge)) if pathway_graphs[pathway].nodes[edge[n]]['TYPE'] == 'Reaction']
            if len(reaction) == 2:
                # reaction 1 --> reaction 2
                for m1 in pathway_graphs[pathway].in_edges(edge[0]):
                    if pathway_graphs[pathway].nodes[m1[0]]['TYPE'] == 'Reaction':
                        continue
                    for m2 in pathway_graphs[pathway].out_edges(edge[1]):
                        if pathway_graphs[pathway].nodes[m2[1]]['TYPE'] == 'Reaction' or m1[0] == m2[1]:
                            continue
                        #G.add_edge(m1[0], m2[1], type=pathway_graphs[pathway].edges[edge]['TYPE'], reaction='')

            if len(reaction) == 1:
                reaction = reaction[0]
                metabolite = [n for n in range(len(edge)) if pathway_graphs[pathway].nodes[edge[n]]['TYPE'] == 'Metabolite'][0]
                if reaction == 0:
                    # reaction --> metabolite; get all metabolites inputted into the reaction (in_edges)
                    for m in pathway_graphs[pathway].in_edges(edge[reaction]):
                        if pathway_graphs[pathway].nodes[m[0]]['TYPE'] == 'Reaction' or m[0] == edge[metabolite]:
                            continue
                        G.add_edge(m[0], edge[metabolite], type=pathway_graphs[pathway].edges[edge]['TYPE'], reaction=pathway_graphs[pathway].nodes[edge[reaction]]['displayName'])
                if reaction == 1:
                    # metabolite --> reaction; get all metabolites outputted from the reaction (out_edges)
                    for m in pathway_graphs[pathway].out_edges(edge[reaction]):
                        if pathway_graphs[pathway].nodes[m[1]]['TYPE'] == 'Reaction' or edge[metabolite] == m[1]:
                            continue
                        G.add_edge(edge[metabolite], m[1], type=pathway_graphs[pathway].edges[edge]['TYPE'], reaction=pathway_graphs[pathway].nodes[edge[reaction]]['displayName'])
        
        nodes, edges = [{'data': G.nodes[node], 'group': 'nodes'} for node in G.nodes], [{'data': {'source': edge[0], 'target': edge[1], 'type': G.edges[edge]['type'], 'reaction': G.edges[edge]['reaction']}, 'group': 'edges'} for edge in G.edges]
        pathway_name = pathway+' ('+str(pathways[1].count(pathway))+')'
        pathway_names.extend([pathway_name])
        pathway_data[pathway_name] = {'nodes': nodes, 'edges': edges}

print(json.dumps({'pathway_data': pathway_data, 'pathway_list': pathway_names}))


# Add in motif and sub-motif options
# Work on saving session, results, SVG or PNG images etc.

# Preview entire plot in small pic in corner like Reactome?
# node popup menu with details
# proteomics, RNAseq
# link pathways together
# missing metabolites in sample cohorts (COLMARq; p-values)
# annotate pathway name over top of pathway sub-graph
#  - select metabolite, highlight pathways in yellow OR clear list and only show the pathways it is in
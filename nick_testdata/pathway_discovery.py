import pandas, numpy, csv, sys, math
import networkx as nx
sys.path.append('../scripts')
from neo4j_connect import ReactomeServer
Reactome = ReactomeServer()

with open('nick_testdata.xlsx', 'rb') as xlsx:
    data = pandas.read_excel(xlsx)

title = 'CF_NTvsBC6'
#print(data.keys())
data_set = numpy.array([data[title], data['Unnamed: 17'], data['Unnamed: 18']])

names_to_colmarid = {}
with open('../updating_Reactome_molecules/name_to_colmarid.csv') as csvin:
    csvin.readline()
    for row in csv.reader(csvin):
        names_to_colmarid[row[0]] = row[1]

metabolites = {}
for i in range(len(data_set[0])):
    if data_set[0][i] == 'METABOLITES':
        continue
    metabolite, p_value, fold_change = data_set[0][i], -math.log(data_set[1][i]), data_set[2][i]
    if metabolite in names_to_colmarid:
        colmar = names_to_colmarid[metabolite]
        metabolites[metabolite] = {'COLMAR': colmar, 'Fold Change': fold_change, 'P-Value': p_value}

# 1. Search for pathways with any of the inputted metabolites
pathways = [{}, [], [], []]
with Reactome._driver.session() as db:
    for metabolite in metabolites:
        pathways[0][metabolite] = set()
        colmar = metabolites[metabolite]['COLMAR']
        fold_change = metabolites[metabolite]['Fold Change']
        for (m, p) in db.run('MATCH (m:ReferenceMolecule)-[me:referenceEntity]-(e:PhysicalEntity)-[er:input]-(r:Reaction)-[rp:hasEvent]-(p:Pathway) WHERE "'+colmar+'" IN m.COLMAR AND p.speciesName = "Homo sapiens" RETURN ([m.displayName, p.displayName])').value():
            pathways[0][metabolite].add(p)
        #print(metabolite, fold_change)

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
    if count > 3:
        pathways[3].extend([p])

# 6. Search for all reactions in each of those pathways
pathway_graphs = {}
with Reactome._driver.session() as db:
    # https://neo4j.com/docs/api/python-driver/current/types/graph.html
    entities_to_metabolites = {}
    #for pathway in pathways[3]:
    for pathway in ["Citric acid cycle (TCA cycle)"]:
        pathway_graphs[pathway] = nx.DiGraph()

        command = 'MATCH (p:Pathway) WHERE p.displayName = "'+pathway+'" AND p.speciesName = "Homo sapiens"'
        command = command+'\n'+'MATCH (R1:Reaction)-[pe:precedingEvent]-(R2:Reaction) WHERE (R1)-[:hasEvent]-(p)  AND (R2)-[:hasEvent]-(p) '
        command = command+'\n'+'RETURN R1, R2, pe'

        reactions = db.run(command).graph()

        if reactions.nodes:
            print(pathway, '\n')
            for node in reactions.nodes:
                if 'Reaction' in list(node.labels): 
                    node_dict = dict(node)
                    node_dict['labels'], node_dict['TYPE'] = list(node.labels), 'Reaction'
                    pathway_graphs[pathway].add_nodes_from([(node['dbId'], node_dict)])

                    command = 'MATCH (p:Pathway) WHERE p.displayName = "'+pathway+'" AND p.speciesName = "Homo sapiens"'
                    command = command+'\n'+'MATCH (m:ReferenceMolecule)-[me:referenceEntity]-(e:PhysicalEntity)-[er]-(r:Reaction)-[rp:hasEvent]-(p) WHERE r.dbId = '+str(node['dbId'])+' RETURN *'
                    reaction = db.run(command).graph()
                    
                    metabolites = [n['displayName'] for n in reaction.nodes if 'ReferenceMolecule' in list(n.labels)]
                    
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

        G = nx.DiGraph()
        for node in pathway_graphs[pathway].nodes:
            G.add_node(node, name=pathway_graphs[pathway].nodes[node]['displayName'], type=pathway_graphs[pathway].nodes[node]['TYPE'])
        for edge in pathway_graphs[pathway].edges:
            G.add_edge(edge[0], edge[1], type=pathway_graphs[pathway].edges[edge]['TYPE'])
        nx.write_graphml(G, 'TCA.graphml')
            

exit()
# 7. Construct new graph that contains pathway reactions as metabolites involved in those reactions that are in order of the pathway order
G = nx.Graph()
for metabolite in pathways[0]:
    G.add_node(metabolite, type='metabolite', fold_change=metabolites[metabolite]['Fold Change'], pvalue=metabolites[metabolite]['P-Value'])
    for pathway in pathways[0][metabolite]:
        if pathway in pathways[3]:
            if pathway not in G.nodes:
                G.add_node(pathway, type='pathway')
            G.add_edge(metabolite, pathway)

nx.write_graphml(G, title+'.graphml')
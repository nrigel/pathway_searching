import pandas, numpy, csv, sys, math
import networkx as nx
sys.path.append('../scripts')
from neo4j_connect import ReactomeServer
Reactome = ReactomeServer()

with open('nick_testdata.xlsx', 'rb') as xlsx:
    data = pandas.read_excel(xlsx)

data_set = numpy.array([data['NT_CFvsWT'], data['Unnamed: 1'], data['Unnamed: 2']])

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


pathways = [{}, [], [], []]
with Reactome._driver.session() as db:
    for metabolite in metabolites:
        pathways[0][metabolite] = set()
        colmar = metabolites[metabolite]['COLMAR']
        fold_change = metabolites[metabolite]['Fold Change']
        for (m, p) in db.run('MATCH (m:ReferenceMolecule)-[me:referenceEntity]-(e:PhysicalEntity)-[er:input]-(r:Reaction)-[rp:hasEvent]-(p:Pathway) WHERE "'+colmar+'" IN m.COLMAR AND p.speciesName = "Homo sapiens" RETURN ([m.displayName, p.displayName])').value():
            pathways[0][metabolite].add(p)
        #print(metabolite, fold_change)

for metabolite in pathways[0]:
    for pathway in pathways[0][metabolite]:
        pathways[1].extend([pathway])

pathways[1] = sorted(pathways[1], key=lambda x: -pathways[1].count(x))

for p in pathways[1]:
    if p not in pathways[2]:
        pathways[2].extend([p])

for p in pathways[2]:
    count = pathways[1].count(p)
    if count > 3:
        pathways[3].extend([p])

G = nx.Graph()
for metabolite in pathways[0]:
    G.add_node(metabolite, type='metabolite', fold_change=metabolites[metabolite]['Fold Change'], pvalue=metabolites[metabolite]['P-Value'])
    for pathway in pathways[0][metabolite]:
        if pathway in pathways[3]:
            if pathway not in G.nodes:
                G.add_node(pathway, type='pathway')
            G.add_edge(metabolite, pathway)

nx.write_graphml(G, 'pathway_results.graphml')
import xml.etree.ElementTree as ET
import networkx as nx

filename = 'PW002044.owl'

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
    xrefs = [G.nodes[e[1]]['resource'] for e in G.out_edges(pathway) if G.nodes[e[1]]['tag'] == 'xref']
    for x in range(len(xrefs)):
        if 'Reference/' in xrefs[x]:
            xrefs[x] = xrefs[x][len('Reference/'):]
    print(xrefs)
# pathwayOrder --> BiochemicalPathwayStep --> stepConversion --> BioChemical Reaction: left, right, name, spontaneous?, stoichiometry, conversionDirection
# make reaction SVG images using OpenBabel, RDKit?
exit()

# Add pathway order to Graph
for parent in TAGs['Pathway']:
    displayName, organism = None, None
    for child in parent:
        if child.get('tag') == 'displayName':
            displayName = child.text
        if child.get('tag') == 'organism':
            organism = child.get('resource')
    print(displayName, organism)
    G.add_nodes_from([(i, elem_data[parent])])
    pathwayOrder = [c for c in parent if elem_data[int(c)]['tag'] == 'pathwayOrder']
    for child in pathwayOrder:
        ID = elem_data[int(child)]['resource'][1:]

print(list(tags.keys()))
        
# Add SmallMoleculeReference to Graph
for elem in tags['SmallMoleculeReference']:
    #print(elem_data[elem])
    
    
    break
    
#print(tags.keys())


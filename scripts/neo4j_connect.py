from neo4j import GraphDatabase
import networkx as nx

# set parameters for Bolt port, user, and password
dbsettings = {'NMRT': {'port': '7687', 'user': 'neo4j', 'password': 'olivia05'},
          'Reactome': {'port': '7688', 'user': 'neo4j', 'password': 'olivia05'}}

# define both DB server classes with test functions
class NMRTServer(object):

    def __init__(self, uri='bolt://127.0.0.1:'+dbsettings['NMRT']['port']+'/', user=dbsettings['NMRT']['user'], password=dbsettings['NMRT']['password']):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
        
    def close(self):
        self._driver.close()

    def test(self):
        with self._driver.session() as session:
            mets = session.run('MATCH (m:metabolites) RETURN count(m)').value()
        print('NMRT loaded with '+str(mets[0])+' metabolites.')
        return 

    def metabolitegraph(self):
        if 'metabolites' not in dir(self):
            with self._driver.session() as db:
                self.metabolites = {}
                self.metabolites['Neo4j'] = db.run('MATCH (m:metabolites)-[r:metaboliteNode]-(M:metaboliteNodes) RETURN *').graph()
                self.metabolites['networkx'] = nx.Graph()
                self.metabolites['spinsystems'] = {}
                for node in self.metabolites['Neo4j'].nodes:
                    if list(node.labels)[0] == 'metaboliteNodes':
                        self.metabolites['networkx'].add_node(node['Metabolite']+'_'+str(node['Node']), neo4j=node)
                        if node['Metabolite'] not in self.metabolites['spinsystems']:
                            self.metabolites['spinsystems'][node['Metabolite']] = set()
                        self.metabolites['spinsystems'][node['Metabolite']].add(node['SpinSystem'])
                    else:
                        self.metabolites['networkx'].add_node(node['Compound'], neo4j=node)
                for relationship in self.metabolites['Neo4j'].relationships:
                    self.metabolites['networkx'].add_edge(relationship.start_node['Compound'], relationship.end_node['Metabolite']+'_'+str(relationship.end_node['Node']), neo4j=relationship)
        return self.metabolites


class ReactomeServer(object):
    
    def __init__(self, uri='bolt://127.0.0.1:'+dbsettings['Reactome']['port']+'/', user=dbsettings['Reactome']['user'], password=dbsettings['Reactome']['password']):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
        
    def close(self):
        self._driver.close()

    def test(self):
        with self._driver.session() as session:
            mets = session.run('MATCH (m:Metabolite) RETURN count(m)').value()
        print('Reactome loaded with '+str(mets[0])+' metabolites.')
        return 

    def generategraphs(self):
        if 'graphs' not in dir(self): # pull metabolite to pathways, metabolite, submotif, and motif grpahs from Neo4j DB
            # https://neo4j.com/docs/api/python-driver/current/types/graph.html
            self.graphs = {'submotifs': {}, 'motifs': {}, 'metabolites': {}}
            with self._driver.session() as db:
                
                self.graphs['metabolites_to_pathways'] = {}
                self.graphs['metabolites_to_pathways']['Neo4j'] = db.run('MATCH (m:Metabolite)-[r:metabolizedIn]-(p) RETURN *').graph()
                self.graphs['metabolites_to_pathways']['networkx'] = nx.Graph()
                for node in self.graphs['metabolites_to_pathways']['Neo4j'].nodes:
                    if 'Metabolite' in list(node.labels):
                        self.graphs['metabolites_to_pathways']['networkx'].add_node(node['SMILES'], labels=node.labels, neo4j=node)
                    else:
                        self.graphs['metabolites_to_pathways']['networkx'].add_node(node['dbId'], labels=node.labels, neo4j=node)
                for relationship in self.graphs['metabolites_to_pathways']['Neo4j'].relationships:
                    self.graphs['metabolites_to_pathways']['networkx'].add_edge(relationship.start_node['SMILES'], relationship.end_node['dbId'], neo4j=relationship)

                for shell in range(1, 5):
                    self.graphs['submotifs'][shell] = {}
                    self.graphs['submotifs'][shell]['Neo4j'] = db.run('MATCH (s:submotif_'+str(shell)+')-[r:submotif'+str(shell)+']-(m:Metabolite) RETURN *').graph()
                    self.graphs['submotifs'][shell]['networkx'] = nx.Graph()
                    for node in self.graphs['submotifs'][shell]['Neo4j'].nodes:
                        self.graphs['submotifs'][shell]['networkx'].add_node(node['SMILES'], labels=node.labels, neo4j=node)
                    for relationship in self.graphs['submotifs'][shell]['Neo4j'].relationships:     
                        self.graphs['submotifs'][shell]['networkx'].add_edge(relationship.start_node['SMILES'], relationship.end_node['SMILES'], metabolite=relationship.start_node['SMILES'], neo4j=relationship)
                
                for shell in range(0, 3):
                    self.graphs['motifs'][shell] = {}
                    self.graphs['motifs'][shell]['Neo4j'] = db.run('MATCH (s:motif_'+str(shell)+')-[r:motif'+str(shell)+']-(m:Metabolite) RETURN *').graph()
                    self.graphs['motifs'][shell]['networkx'] = nx.Graph()
                    for node in self.graphs['motifs'][shell]['Neo4j'].nodes:
                        self.graphs['motifs'][shell]['networkx'].add_node(node['SMILES'], labels=node.labels, neo4j=node)
                    for relationship in self.graphs['motifs'][shell]['Neo4j'].relationships:     
                        self.graphs['motifs'][shell]['networkx'].add_edge(relationship.start_node['SMILES'], relationship.end_node['SMILES'], metabolite=relationship.start_node['SMILES'], neo4j=relationship)

    def pathwaymatch(self, SMILES, structuretype, shell):
        self.generategraphs()
        P = self.graphs['metabolites_to_pathways']['networkx']
        results = {}
        if structuretype == 'motifs' or structuretype == 'submotifs':
            G = self.graphs[structuretype][shell]['networkx']
            if SMILES in G.nodes:
                for edge in G.edges(SMILES):
                    for p in P.edges(G.edges[edge]['metabolite']):
                        pathway = [n for n in p if n != G.edges[edge]['metabolite']][0]
                        name, species = P.nodes[pathway]['neo4j']['displayName'], P.nodes[pathway]['neo4j']['speciesName']
                        if species not in results:
                            results[species] = []
                        results[species].extend([name])
        
        elif structuretype == 'metabolites':
            for p in P.edges(SMILES):
                pathway = [n for n in p if n != SMILES][0]
                name, species = P.nodes[pathway]['neo4j']['displayName'], P.nodes[pathway]['neo4j']['speciesName']
                if species not in results:
                    results[species] = []
                results[species].extend([name])

        return results

    def pathwayvis(self):
        # here, we could edit the returned graphs to produce a visualization
        pass
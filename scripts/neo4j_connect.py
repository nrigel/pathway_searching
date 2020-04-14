from neo4j import GraphDatabase

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

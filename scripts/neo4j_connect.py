from neo4j import GraphDatabase

dbsettings = {'NMRT': {'port': '7687', 'user': 'neo4j', 'password': 'olivia05'},
          'Reactome': {'port': '7688', 'user': 'neo4j', 'password': 'olivia05'}}


class NMRTServer(object):

    def __init__(self, uri='bolt://127.0.0.1:'+dbsettings['NMRT']['port']+'/', user=dbsettings['NMRT']['user'], password=dbsettings['NMRT']['password']):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))
    
        print(self.test())

    def close(self):
        self._driver.close()

    def test(self):
        with self._driver.session() as session:
            return session.run('MATCH (m:metabolites) RETURN count(m)').value()


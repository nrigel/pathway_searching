dbsettings = {'NMRT': {'port': '7687', 'user': 'neo4j', 'password': 'olivia05'},
      'PathBank': {'port': '7690', 'user': 'neo4j', 'password': 'olivia05'}}
from neo4j import GraphDatabase

class Neo4jServer(object):
    def __init__(self, Server):
        uri = 'bolt://127.0.0.1:'+dbsettings[Server]['port']+'/'
        user = dbsettings[Server]['user']
        password = dbsettings[Server]['password']
        self._driver = GraphDatabase.driver(uri, auth=(user, password))   
    def close(self):
        self._driver.close()
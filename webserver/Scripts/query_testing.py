import time

from neo4j_connect import Neo4jServer

NMRT, PathBank = Neo4jServer('NMRT'), Neo4jServer('PathBank')

with PathBank._driver.session() as db:
    start = time.time()
    db.run('MATCH (p:Pathway {SMPDB_ID: "SMP0000067"}) RETURN p').value()
    print(int(round(time.time()-start, 3)*1000), 'ms to match pathway only')

    start = time.time()
    db.run('MATCH (p:Pathway {SMPDB_ID: "SMP0000067"})-[he:hasEvent]-(r)-[er:input|output]-(m) RETURN ([apoc.map.removeKey(m {.*}, "SMPDB_ID"), type(er), r, he, p])').value()
    print(int(round(time.time()-start, 3)*1000), 'ms to match everything')

    start = time.time()
    result = db.run('MATCH (p:Pathway {SMPDB_ID: "SMP0000067"})-[he:hasEvent]-(r)-[er:input|output]-(m) RETURN ([[x in keys(m) WHERE not x in ["SMPDB_ID"]| [x, m[x] ] ], type(er), r, he, p])').value()
    print(int(round(time.time()-start, 3)*1000), 'ms to match everything')

    start = time.time()
    result = db.run('MATCH (p:Pathway {SMPDB_ID: "SMP0000067"})-[he:hasEvent]-(r)-[er:input|output]-(m) RETURN ([ [x in keys(m) WHERE not x in ["SMPDB_ID"]| [x, m[x] ] ], type(er), [x in keys(r) | [x, r[x] ] ], he, [x in keys(p) | [x, p[x] ] ]])').value()
    print(int(round(time.time()-start, 3)*1000), 'ms to match everything')

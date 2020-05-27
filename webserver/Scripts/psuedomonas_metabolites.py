from neo4j_connect import Neo4jServer

NMRT, PathBank = Neo4jServer('NMRT'), Neo4jServer('PathBank')

metabolites = set()

with PathBank._driver.session() as db:
    pathways = db.run('MATCH (p:Pathway) WHERE p.Species = "Pseudomonas aeruginosa" RETURN p.SMPDB_ID').value()
    for pathway in pathways:
        for m in db.run('MATCH (m:Metabolite) WHERE "'+pathway+'" IN m.SMPDB_ID RETURN ([m.PW_ID, m.HMDB_ID])').value():
            metabolites.add((m[0], m[1]))

import csv

with open('Pseudomonas_metabolites.csv', 'w') as csvout:
    writer = csv.writer(csvout)
    writer.writerow(['PW_ID', 'HMDB_ID'])
    for m in metabolites:
        writer.writerow([m[0], m[1]])
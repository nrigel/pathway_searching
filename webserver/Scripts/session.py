#!/Users/nick/Documents/GitHub/motif_builder/py35env/bin/python

print("Content-type:text/html\r\n\r\n") # print HTML header

from neo4j_connect import Neo4jServer
PathBank = Neo4jServer('PathBank')

import cgi, cgitb, json, traceback

from session_save import SaveSession

form = cgi.FieldStorage() # Create instance of FieldStorage 
saveinfo = json.loads(form.getvalue('saveinfo'))
cached_results = json.loads(form.getvalue('cached_results'))

try:
    with PathBank._driver.session() as db:
        if saveinfo['save']:
            saveinfo = SaveSession(db, cached_results, saveinfo)
        
        if saveinfo['load']:
            pass

    PathBank.close()

    print(json.dumps(cached_results))

except Exception as e:
    print(traceback.format_exc())
    exit()
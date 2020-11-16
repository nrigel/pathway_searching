#!/Users/nick/Documents/GitHub/motif_builder/py35env/bin/python

try:
    from neo4j_connect import Neo4jServer
    PathBank = Neo4jServer('PathBank')

    import sys, json, traceback

    from session_save import SaveSession
    

    form = json.loads(sys.argv[1]) # load the form data from php
    saveinfo = json.loads(form['saveinfo'])
    cached_results = json.loads(form['cached_results'])

    with PathBank._driver.session() as db:
        if 'save' in saveinfo and saveinfo['save']:
            saveinfo = SaveSession(db, cached_results, saveinfo)
            saveinfo['save'] = False
        
        if 'load' in saveinfo and saveinfo['load']:
            if 'session_id' in saveinfo and saveinfo['session_id']:
                session = db.run('MATCH (n:Sessions {Name: "'+saveinfo['session_id']+'"}) RETURN n').value()
                if session:
                    session = dict(session[0])
                    saveinfo['saveStateOrder'] = json.loads(session['saveStateOrder'])
                    saveinfo['saveStateNames'] = json.loads(session['saveStateNames'])
                    saveinfo['saveStates'] = json.loads(session['saveStates'])
                    saveinfo['lastSave'] = json.loads(session['lastSave'])
                    sp = saveinfo['saveStateOrder'][0]
                    if saveinfo['savestate'] and saveinfo['savestate'] in saveinfo['saveStateNames']:
                        sp = saveinfo['saveStateNames'][saveinfo['savestate']]
                    cached_results = json.loads(session[sp])
                    saveinfo['savestate'] = saveinfo['saveStates'][sp]

                    # add SVGs
                    from svg_drawer import Draw as SVGDraw
                    for opt in cached_results["pathway_data"]:
                        for pathway in cached_results["pathway_data"][opt]:
                            for node in cached_results["pathway_data"][opt][pathway]['nodes']:
                                if 'SMILES' not in node['data']:
                                    continue
                                smiles = node['data']['SMILES']
                                if '_' in smiles:
                                    smiles = smiles.split('_')[0]
                                node['data']['SVG'] = SVGDraw(smiles)
                    
            saveinfo['load'] = False

    PathBank.close()

    cached_results['saveinfo'] = saveinfo
    print(json.dumps(cached_results))

except Exception as e:
    print(traceback.format_exc())
    exit()
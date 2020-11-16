#!/Users/nick/Documents/GitHub/motif_builder/py35env/bin/python

import json

def SaveSession(db, session_data, saveinfo):
    
    def getCurrentTime():
        import datetime

        d, t = str(datetime.datetime.now()).split()
        d, t = '/'.join([str(int(i)) for i in d.split('-')[1:]]+[d.split('-')[0]]), ':'.join(t.split(':')[:2])
        return 'Save point '+t+' '+d

    if not saveinfo['savestate']: # set savestate name
        saveinfo['savestate'] = getCurrentTime()

    if not saveinfo['session_id']: # create session node
        from random import choice
        from string import ascii_uppercase
        
        used_names = db.run('MATCH (n:Sessions) RETURN n.Name').value()
        saveinfo['session_id'] = str(len(used_names)+1).zfill(4)+'-'+''.join(choice(ascii_uppercase) for i in range(10))
        db.run('CREATE (n:Sessions {Name: "'+saveinfo['session_id']+'", saveStates: "{}", saveStateOrder: "[]", lastSave: "{}"})') 

    saveinfo['saveStates'], saveinfo['saveStateOrder'], saveinfo['lastSave'] = db.run('MATCH (n:Sessions {Name: "'+saveinfo['session_id']+'"}) RETURN ([n.saveStates, n.saveStateOrder, n.lastSave])').value()[0]
    saveinfo['saveStates'], saveinfo['saveStateOrder'], saveinfo['lastSave'] = json.loads(saveinfo['saveStates']), json.loads(saveinfo['saveStateOrder']), json.loads(saveinfo['lastSave'])
    saveinfo['saveStateNames'] = {saveinfo['saveStates'][key]: key for key in saveinfo['saveStates']}

    if saveinfo['savestate'] not in saveinfo['saveStateNames']:
        saveinfo['saveStateNames'][saveinfo['savestate']] = 'sp'+str(len(saveinfo['saveStateNames']))
        saveinfo['saveStates'][saveinfo['saveStateNames'][saveinfo['savestate']]] = saveinfo['savestate']

    if saveinfo['saveStateNames'][saveinfo['savestate']] in saveinfo['saveStateOrder']:
        saveinfo['saveStateOrder'].remove(saveinfo['saveStateNames'][saveinfo['savestate']])

    saveinfo['saveStateOrder'].insert(0, saveinfo['saveStateNames'][saveinfo['savestate']])
    saveinfo['lastSave'][saveinfo['savestate']] = getCurrentTime()
    
    # clean up our session data: remove any SVGs
    for opt in session_data["pathway_data"]:
        for pathway in session_data["pathway_data"][opt]:
            for node in session_data["pathway_data"][opt][pathway]['nodes']:
                node['data']['SVG'] = None

    db.run('MATCH (n:Sessions {Name: "'+saveinfo['session_id']+'"}) SET n.'+saveinfo['saveStateNames'][saveinfo['savestate']]+' = '+json.dumps(json.dumps(session_data))+', n.saveStates = '+json.dumps(json.dumps(saveinfo['saveStates']))+', n.saveStateOrder = '+json.dumps(json.dumps(saveinfo['saveStateOrder']))+', n.saveStateNames = '+json.dumps(json.dumps(saveinfo['saveStateNames']))+', n.lastSave = '+json.dumps(json.dumps(saveinfo['lastSave'])))

    return saveinfo
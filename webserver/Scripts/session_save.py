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
        db.run('CREATE (n:Sessions {Name: "'+saveinfo['session_id']+'", saveStates: "{}", saveStateOrder: "[]"})') 

    saveStates, saveStateOrder = db.run('MATCH (n:Sessions {Name: "'+saveinfo['session_id']+'"}) RETURN ([n.saveStates, n.saveStateOrder])').value()[0]
    saveStates, saveStateOrder = json.loads(saveStates), json.loads(saveStateOrder)
    saveinfo['saveStateNames'] = {saveStates[key]: key for key in saveStates}

    if saveinfo['savestate'] not in saveinfo['saveStateNames']:
        saveinfo['saveStateNames'][saveinfo['savestate']] = 'sp'+str(len(saveinfo['saveStateNames']))
        saveStates[saveinfo['saveStateNames'][saveinfo['savestate']]] = saveinfo['savestate']

    if saveinfo['saveStateNames'][saveinfo['savestate']] in saveStateOrder:
        saveStateOrder.remove(saveinfo['saveStateNames'][saveinfo['savestate']])

    saveStateOrder.insert(0, saveinfo['saveStateNames'][saveinfo['savestate']])

    # clean up our session data: remove any SVG
    print(json.dumps(session_data))

    db.run('MATCH (n:Sessions {Name: "'+saveinfo['session_id']+'"}) SET n.'+saveinfo['saveStateNames'][saveinfo['savestate']]+' = '+json.dumps(json.dumps(session_data))+', n.saveStates = '+json.dumps(json.dumps(saveStates))+', n.saveStateOrder = '+json.dumps(json.dumps(saveStateOrder)))

    return saveinfo
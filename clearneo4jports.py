# automatically find old processes running on the Neo4j Bolt ports and kill them

import os

pids = set()
for port in [7687, 7688, 7476, 7689]:

    response = os.popen('lsof -i :'+str(port)).read()

    if response:
        for process in response.split('\n')[1:]:
            if process.split():
                pids.add(process.split()[1])
if pids:
    os.popen('kill '+' '.join(pids))
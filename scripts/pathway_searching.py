
# initialize Neo4j Databases
from neo4j_connect import NMRTServer, ReactomeServer
NMRT, Reactome = NMRTServer(), ReactomeServer()

class PathwayMapper(object):

    def __init__(self, matching_cutoff=5, replicates=3):
        # store our initial variables
        self.matching_cutoff, self.replicates = matching_cutoff, replicates

        self.getcsaverages()

    def getcsaverages(self):
        # collect average standard deviations for C and H chemical shifts following the replicate cutoff
        a = ['null as Cavg'+str(i)+', null as Havg'+str(i) for i in range(1, 5)]
        a = 'WITH '+', '.join(a)
        b = ['MATCH (n:submotif_'+str(i)+') WHERE toInteger(split(n.CS, " ")[0]) >= '+str(self.replicates)+' WITH sum(toFloat(split(n.CS, " ")[2]))/count(split(n.CS, " ")) AS Cavg'+str(i)+', sum(toFloat(split(n.CS, " ")[4]))/count(split(n.CS, " ")) AS Havg'+str(i) for i in range(1, 5)]
        
        
#RETURN Cavg1, Havg1


    def spinsystem(self, ss):
        # first search for sub-motif matches
        for node in ss:
            for C in ss[node]:
                self.submotifmatcher(C, ss[node][C])
            break

    def submotifmatcher(self, C, H):
        with NMRT._driver.session() as db:
            # submotif.CS = [replicates, Cavg, Cstd, Havg, Hstd]
            print(db.run('MATCH (s:submotif_1) RETURN s.CS').value())


# example spin-systems
arginine_ss = {8: {43.2: 3.2364}, 4: {57.003: 3.7638}, 7: {26.576: 1.6795}, 6: {30.281: 1.909}}

pm = PathwayMapper()
pm.spinsystem(arginine_ss)
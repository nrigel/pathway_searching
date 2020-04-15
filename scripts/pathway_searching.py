
# initialize Neo4j Databases
from neo4j_connect import NMRTServer, ReactomeServer
NMRT, Reactome = NMRTServer(), ReactomeServer()

class PathwayMapper(object):

    def __init__(self, matching_cutoff=5, replicates=3):
        # store our initial variables
        self.matching_cutoff, self.replicates = matching_cutoff, replicates
        self.r = (1, 5)
        self.csaverages = self.getcsaverages()

    def getcsaverages(self):
        # collect average standard deviations for C and H chemical shifts following the replicate cutoff
        a = 'WITH '+', '.join(['null as Cavg'+str(i)+', null as Havg'+str(i) for i in range(self.r[0], self.r[1])])
        b = ['Cavg'+str(i) for i in range(self.r[0], self.r[1])]+['Havg'+str(i) for i in range(self.r[0], self.r[1])]
        c = []
        for i in range(self.r[0], self.r[1]):
            C = '\nMATCH (n:submotif_'+str(i)+') WHERE toInteger(split(n.CS, " ")[0]) >= '+str(self.replicates)+' WITH sum(toFloat(split(n.CS, " ")[2]))/count(split(n.CS, " ")) AS Cavg'+str(i)+', sum(toFloat(split(n.CS, " ")[4]))/count(split(n.CS, " ")) AS Havg'+str(i)
            for B in b:
                if B != 'Cavg'+str(i) and B != 'Havg'+str(i):
                    C = C+', '+B
            c.extend([C])
        d = '\nRETURN '+', '.join(['Cavg'+str(i)+', Havg'+str(i) for i in range(self.r[0], self.r[1])])
        with NMRT._driver.session() as db:
            avgs = db.run(a+''.join(c)+d).values()[0]
        AVGs = {}
        for i in range(0, len(avgs), 2):
            AVGs[len(AVGs)+1] = {avgs[i]: avgs[i+1]}
        return AVGs

    def submotifmatcher(self, C, H):
        with NMRT._driver.session() as db:
            # submotif.CS = [replicates, Cavg, Cstd, Havg, Hstd]
            matches = {}
            for i in range(self.r[0], self.r[1]):
                Cavg = list(self.csaverages[i].keys())[0]
                Havg = self.csaverages[i][Cavg]
                Cavg, Havg = 1.9005, 0.1165
                a = '\n'.join(['WITH '+str(C)+' AS Cexp, '+str(H)+' AS Hexp, '+str(Cavg)+' AS Cavg, '+str(Havg)+' AS Havg, '+str(self.matching_cutoff)+' as cutoff',
                                'MATCH (s:submotif_'+str(i)+') RETURN', 'CASE', 'WHEN toInteger(split(s.CS, " ")[0]) >= '+str(self.replicates),
                                'THEN CASE', 'WHEN toInteger(split(s.CS, " ")[2]) > 0 AND toInteger(split(s.CS, " ")[4]) > 0',
                                'THEN CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2 <= cutoff',
                                'THEN [s.SMILES, (abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2] END',
                                'WHEN toInteger(split(s.CS, " ")[2]) = 0 AND toInteger(split(s.CS, " ")[4]) > 0',
                                'THEN CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2 <= cutoff',
                                'THEN [s.SMILES, (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2] END',
                                'WHEN toInteger(split(s.CS, " ")[2]) > 0 AND toInteger(split(s.CS, " ")[4]) = 0', 'THEN CASE',
                                'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/Havg)/2 <= cutoff',
                                'THEN [s.SMILES, (abs(Cexp-toFloat(split(s.CS, " ")[2]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-Havg)/Havg)/2] END',
                                'END', 'ELSE CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/Havg)/2 <= cutoff',
                                'THEN [s.SMILES, (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-Havg)/toFloat(split(s.CS, " ")[3]))/2] END END AS result'])
                matches[i] = {m[0][0]: m[0][1] for m in db.run(a).values() if m[0] != None}
            return matches

    def spinsystem(self, ss):
        # first search for sub-motif matches
        print(self.csaverages)
        for node in ss:
            for C in ss[node]:
                matches = self.submotifmatcher(C, ss[node][C])
                
                print(node, len(matches[2]))

    
# example spin-systems
arginine_ss = {8: {43.2: 3.2364}, 4: {57.003: 3.7638}, 7: {26.576: 1.6795}, 6: {30.281: 1.909}}

pm = PathwayMapper()
pm.spinsystem(arginine_ss)

NMRT.close()
Reactome.close()
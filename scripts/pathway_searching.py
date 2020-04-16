# initialize Neo4j Databases
from neo4j_connect import NMRTServer, ReactomeServer
NMRT, Reactome = NMRTServer(), ReactomeServer()

# import libraries for bipartite z-scoring
from munkres import Munkres
import numpy as np
import math

import networkx as nx

class PathwayMapper(object):

    def __init__(self, matching_cutoff=5, replicates=3):
        # store our initial variables
        self.matching_cutoff, self.replicates = matching_cutoff, replicates
        self.r = (1, 5)
        self.csaverages = self.getcsaverages()

    def getcsaverages(self):
        # collect average standard deviations for C and H chemical shifts following the replicate cutoff
        AVGs = {}
        ranges = {'submotif': self.r, 'motifnode': (0, 3)}
        with NMRT._driver.session() as db:
            for m in ranges:
                AVGs[m] = {}
                r = ranges[m]
                a = 'WITH '+', '.join(['null as Cavg'+str(i)+', null as Havg'+str(i) for i in range(r[0], r[1])])
                b = ['Cavg'+str(i) for i in range(r[0], r[1])]+['Havg'+str(i) for i in range(r[0], r[1])]
                c = []
                for i in range(r[0], r[1]):
                    C = '\nMATCH (n:'+m+'_'+str(i)+') WHERE toInteger(split(n.CS, " ")[0]) >= '+str(self.replicates)+' WITH sum(toFloat(split(n.CS, " ")[2]))/count(split(n.CS, " ")) AS Cavg'+str(i)+', sum(toFloat(split(n.CS, " ")[4]))/count(split(n.CS, " ")) AS Havg'+str(i)
                    for B in b:
                        if B != 'Cavg'+str(i) and B != 'Havg'+str(i):
                            C = C+', '+B
                    c.extend([C])                
                d = '\nRETURN '+', '.join(['Cavg'+str(i)+', Havg'+str(i) for i in range(r[0], r[1])])
                avgs = db.run(a+''.join(c)+d).values()[0]
                for i in range(0, len(avgs), 2):
                    AVGs[m][len(AVGs[m])+r[0]] = {avgs[i]: avgs[i+1]}
        return AVGs
    
    def bipartitescoring(self, csdict1, csdict2, method):
        # csdict1 (exp): {node: {C: H}}
        # csdict2 (db): {'CS': {node: {(Cavg, Cstd): (Havg, Hstd)}}, 'replicates': {node: replicates}}
        i, matrix = list(csdict1.keys()), []
        if method == 'zscore':
            j, r = list(csdict2['CS'].keys()), float(sum([csdict2['replicates'][n] for n in csdict2['replicates']]))
        if method == 'rmsd':
            j = list(csdict2.keys())

        for n1 in i:
            temp = []
            C1 = list(csdict1[n1].keys())[0]
            H1 = csdict1[n1][C1]
            for n2 in j:
                if method == 'zscore':
                    C2 = list(csdict2['CS'][n2].keys())[0]
                    H2 = csdict2['CS'][n2][C2]
                    R2 = csdict2['replicates'][n2]
                    e = int((abs(C1-C2[0])/C2[1]+abs(H1-H2[0])/H2[1])/2.0*1000)*(1-R2/r) #scaling factor incorporates weighted average, bipartite matching minimizes cost so we want to multiply by the inverse weight
                if method == 'rmsd':
                    C2 = list(csdict2[n2].keys())[0]
                    H2 = csdict2[n2][C2]
                    e = int(math.sqrt((C1-C2)**2/2.0+((H1-H2)**2)*50.0)*1000)
                temp.append(e) # weight for weighted bipartite matching
            matrix.append(temp)
        m = Munkres()
        indexes = np.array(m.compute(matrix))
        score, assignments = 0, {}
        for ind in indexes:
            n1, n2 = i[ind[0]], j[ind[1]]
            assignments[n1] = n2
            C1 = list(csdict1[n1].keys())[0]
            H1 = csdict1[n1][C1]
            if method == 'zscore':
                C2 = list(csdict2['CS'][n2].keys())[0]
                H2 = csdict2['CS'][n2][C2]
                score += (abs(C1-C2[0])/C2[1]+abs(H1-H2[0])/H2[1])/2.0*csdict2['replicates'][n2]
            if method == 'rmsd':
                C2 = list(csdict2[n2].keys())[0]
                H2 = csdict2[n2][C2]
                score += ((C1-C2)**2+((H1-H2)*10)**2)
        
        if method == 'zscore':
            score = score/r
        if method == 'rmsd':
            score = math.sqrt(score/float(2*len(csdict1)))

        return score, assignments

    def submotifmatcher(self, C, H):
        with NMRT._driver.session() as db:
            # submotif.CS = [replicates, Cavg, Cstd, Havg, Hstd]
            matches = {}
            for i in range(self.r[0], self.r[1]):
                Cavg = list(self.csaverages['submotif'][i].keys())[0]
                Havg = self.csaverages['submotif'][i][Cavg]
                a = '\n'.join(['WITH '+str(C)+' AS Cexp, '+str(H)+' AS Hexp, '+str(Cavg)+' AS Cavg, '+str(Havg)+' AS Havg, '+str(self.matching_cutoff)+' as cutoff',
                                'MATCH (s:submotif_'+str(i)+') RETURN', 'CASE', 'WHEN toInteger(split(s.CS, " ")[0]) >= '+str(self.replicates),
                                'THEN CASE', 'WHEN toInteger(split(s.CS, " ")[2]) > 0 AND toInteger(split(s.CS, " ")[4]) > 0',
                                'THEN CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2 <= cutoff',
                                'THEN [s.SMILES, [(abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2, toInteger(split(s.CS, " ")[0])]] END',
                                'WHEN toInteger(split(s.CS, " ")[2]) = 0 AND toInteger(split(s.CS, " ")[4]) > 0',
                                'THEN CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2 <= cutoff',
                                'THEN [s.SMILES, [(abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/toFloat(split(s.CS, " ")[4]))/2, toInteger(split(s.CS, " ")[0])]] END',
                                'WHEN toInteger(split(s.CS, " ")[2]) > 0 AND toInteger(split(s.CS, " ")[4]) = 0', 'THEN CASE',
                                'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-toFloat(split(s.CS, " ")[3]))/Havg)/2 <= cutoff',
                                'THEN [s.SMILES, [(abs(Cexp-toFloat(split(s.CS, " ")[2]))/toFloat(split(s.CS, " ")[2])+abs(Hexp-Havg)/Havg)/2, toInteger(split(s.CS, " ")[0])]] END',
                                'END', 'ELSE CASE', 'WHEN (abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-toFloat(split(s.CS, " ")[3]))/Havg)/2 <= cutoff',
                                'THEN [s.SMILES, [(abs(Cexp-toFloat(split(s.CS, " ")[1]))/Cavg+abs(Hexp-Havg)/toFloat(split(s.CS, " ")[3]))/2, toInteger(split(s.CS, " ")[0])]] END END AS result'])
                matches[i] = {m[0][0]: m[0][1] for m in db.run(a).values() if m[0] != None}
            return matches

    def motifmatcher(self, ss):
        # pull motifs with same number of nodes
        matches = {}
        with NMRT._driver.session() as db:
            for i in range(0, 3):
                matches[i] = {}
                a = 'MATCH (m:motif_'+str(i)+')--(n:motifnode_'+str(i)+') WHERE size(split(split(m.SMILES, "_")[1], " ")) = '+str(len(ss))+' RETURN ([m.SMILES, n.SMILES, n.CS])'
                motifs = {}
                for result in db.run(a).value():
                    if result[0] not in motifs:
                        motifs[result[0]] = {}
                    motifs[result[0]][result[1]] = result[2]
                for motif in motifs:
                    csdict = {'CS': {}, 'replicates': {}}
                    for node in motifs[motif]:
                        if motifs[motif][node] is None:
                            break
                        r, C, c, H, h = motifs[motif][node].split()
                        if int(r) < self.replicates or float(c) == 0:
                            c = list(self.csaverages['motifnode'][i].keys())[0]
                        if int(r) < self.replicates or float(h) == 0:
                            h = self.csaverages['motifnode'][i][list(self.csaverages['motifnode'][i].keys())[0]]
                        csdict['CS'][int(node.split('_')[1])] = {(float(C), float(c)): (float(H), float(h))}
                        csdict['replicates'][int(node.split('_')[1])] = int(r)
                    zscore, assignments = self.bipartitescoring(ss, csdict, 'zscore')
                    if zscore <= self.matching_cutoff:
                        matches[i][motif] = {'zscore': zscore, 'assignments': assignments, 'CS': csdict}
        return matches
         
    def metabolitematcher(self, ss):
        matches = {}
        metabolites = NMRT.metabolitegraph()
        for metabolite in metabolites['spinsystems']:
            for spinsys in metabolites['spinsystems'][metabolite]:
                if len(spinsys.split()) == len(ss):
                    csdict = {}
                    for node in spinsys.split():
                        neo = metabolites['networkx'].nodes[metabolite+'_'+node]['neo4j']
                        if 'CS' in neo.keys():
                            csdict[int(node)] = {float(neo['CS'].split()[0]): float(neo['CS'].split()[1])}
                    if len(csdict) == len(ss):
                        rmsd, assignments = self.bipartitescoring(ss, csdict, 'rmsd')
                        if rmsd <= self.matching_cutoff:
                            matches[metabolite] = {'rmsd': rmsd, 'assignments': assignments, 'CS': csdict, 'SMILES': metabolites['networkx'].nodes[metabolite]['neo4j']['SMILES']}
        return matches        
            
    def reactomequery(self, submotifmatches, motifmatches, metabolitematches):
        pathway_matches = {'submotifs': {}, 'motifs': {}, 'metabolites': {}}

        Reactome.generategraphs() # generate a few networkx subgraphs
        
        # return pathway match list for each node/motif
        for metabolite in metabolitematches:
            pathway_matches['metabolites'][metabolite] = Reactome.pathwaymatch(metabolitematches[metabolite]['SMILES'], structuretype='metabolites', shell=None)

        for node in submotifmatches:
            pathway_matches['submotifs'][node] = {}
            for shell in submotifmatches[node]:
                pathway_matches['submotifs'][node][shell] = {}
                for submotif in submotifmatches[node][shell]:
                    pathway_matches['submotifs'][node][shell][submotif] = Reactome.pathwaymatch(submotif, structuretype='submotifs', shell=shell)

        for shell in motifmatches:
            pathway_matches['motifs'][shell] = {}
            for motif in motifmatches[shell]:
                # Reactome does not include spin-systems in the SMILES
                pathway_matches['motifs'][shell][motif] = Reactome.pathwaymatch(motif.split('_')[0], structuretype='motifs', shell=shell)
                # motif: {species: pathway list}

        return pathway_matches

    def spinsystem(self, ss):
        metabolite_matches = self.metabolitematcher(ss)
        submotif_matches = {}
        for node in ss:
            for C in ss[node]:
                submotif_matches[node] = self.submotifmatcher(C, ss[node][C])
        motif_matches = self.motifmatcher(ss)    
        pathway_matches = self.reactomequery(submotif_matches, motif_matches, metabolite_matches)
        



        # next steps
            # integrate motif builder/write smaller script than we have and call it from here
            # update reactome database?
                # review what we did last time
                # make new one with all possible compounds
                # then make sub-motif and motif portions that we can search for
                # should be able to build database without chemical shifts
            # write reactome function 
            # write short web-server to demonstrate what we can do
                # include example amino acids with maybe even spin-system graphs
                # drop down menu: select example compound, match for sub-motifs, motifs and search those against Reactome

    
# example spin-systems
arginine_ss = {8: {43.2: 3.2364}, 4: {57.003: 3.7638}, 7: {26.576: 1.6795}, 6: {30.281: 1.909}}

pm = PathwayMapper()
pm.spinsystem(arginine_ss)

NMRT.close()
Reactome.close()
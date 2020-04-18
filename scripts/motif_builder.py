from neo4j_connect import NMRTServer
from matching_functions import MatchingFunctions
import itertools, copy, time

class MotifBuilder(object):

    def __init__(self, matching_cutoff=5, replicates=3):
        self.matching_cutoff, self.replicates = matching_cutoff, replicates
        self.NMRT = NMRTServer()
        self.MatchingFunctions = MatchingFunctions(self.NMRT, self.matching_cutoff, self.replicates)

    def motifBuilderQuery(self, matchingDict):
        # set the motif builder session ID
        session_id, catName = 'testing', 'COLMAR'
        with self.NMRT._driver.session() as session:  
            # delete old session by this name
            session.run('MATCH (n:MBQ) WHERE n.Session_ID = "'+session_id+'" DETACH DELETE n RETURN (n) UNION MATCH (n:MBQNodes) WHERE n.Session_ID = "'+session_id+'" DETACH DELETE n RETURN (n)')

            #command = ['MATCH (s:submotif_2) WHERE (s)<-[:includedIn]-(:Database {Name:"COLMAR"}) AND size(split(s.SMILES, "_")) > 1']
            #command.extend(['RETURN ([avg(toFloat(split(s.CS, " ")[2])), avg(toFloat(split(s.CS, " ")[4]))])'])
            #avgs = session.run('\n'.join(command))

            # get highest cache level -> "N"
            #N_max = max(session.run('MATCH (p:motifProducts) WHERE p.Complete = "True" RETURN (p.N)').value())
            N_max = 5
            if N_max >= len(matchingDict):
                N_ss = len(matchingDict)
                search_for_completemotifs = True
            else:
                N_ss = N_max
                search_for_completemotifs = False
        
            dblabels = session.run('call db.labels()').value()
            command = []
            for label in dblabels:
                rmlabels = ['matchedSubMotif_'+str(session_id)]+['matchedMotifProduct_'+str(i)+'_'+str(session_id) for i in range(2, N_ss)]
                for l in rmlabels:
                    if label[:len(l)] == l:
                        command.extend(['MATCH (n:'+label+') REMOVE n:'+label+' RETURN (n)'])
                        command.extend(['UNION'])
            if len(command) > 0:
                session.run('\n'.join(command[:-1]))

            nodelist = list(matchingDict.keys())
            # set up session
            start = time.time()
            command = ['CREATE (Q: MBQ {Session_ID: "'+str(session_id)+'"})'] # create our session
            # add in cs nodes and matching data
            for node in nodelist:
                #print(len(matchingDict[node]))
                command.extend(['CREATE (n'+str(node)+':MBQNodes {Node:'+str(node)+', Session_ID:"'+str(session_id)+'"})'])
                command.extend(['CREATE (Q)-[:mbqNode]->(n'+str(node)+')'])
            command.extend(['WITH Q'])
            for node in nodelist:
                command.extend(['MATCH (n:MBQNodes) WHERE n.Node = '+str(node)+' AND n.Session_ID = "'+str(session_id)+'"'])
                command.extend(['MATCH (s:submotif_2) WHERE s.SMILES IN ['+', '.join(['"'+s+'"' for s in matchingDict[node].keys()])+']'])
                command.extend(['CREATE (n)-[:matchedSubMotif]->(s)'])
                command.extend(['SET s:matchedSubMotif_'+str(session_id)+'_n'+str(node)+', s:matchedSubMotif_'+str(session_id)])
                command.extend(['RETURN (n)'])
                command.extend(['UNION'])
            session.run('\n'.join(command[:-1]))
            
            motifs = {}
            
            complete = '"False"'
            if 2 == len(nodelist):
                complete = '"True"'
            p1 = session.run('MATCH (s2:matchedSubMotif_'+str(session_id)+')-[r2:motifProductClique]->(p:motifProducts_2)<-[r1:motifProductClique]-(s1:matchedSubMotif_'+str(session_id)+') WHERE p.Complete = '+complete+' AND r1.clique = r2.clique AND ID(r1) <> ID(r2) RETURN ([p.SMILES, r1.clique, s1.SMILES, s2.SMILES, labels(s1), labels(s2)])').value()
            P1 = {}
            for product in p1:
                p, c, s1, s2, l1, l2 = product
                n1 = [int(l.split('_')[-1][1:]) for l in l1 if l[:len('matchedSubMotif_'+str(session_id))] == 'matchedSubMotif_'+str(session_id) and len(l) > len('matchedSubMotif_'+str(session_id))]
                n2 = [int(l.split('_')[-1][1:]) for l in l2 if l[:len('matchedSubMotif_'+str(session_id))] == 'matchedSubMotif_'+str(session_id) and len(l) > len('matchedSubMotif_'+str(session_id))]
                for n in n1:
                    for N in n2:
                        if n == N:
                            continue
                        if p not in P1:
                            P1[p] = []
                        P1[p].extend([tuple(sorted([(n, s1), (N, s2)], key=lambda x: x[0]))])
            
            for i in range(2, N_ss):
                m1 = '['+', '.join(['"'+p+'"' for p in P1.keys()])+']'
                session.run('MATCH (m:motifProducts_'+str(i)+') WHERE m.SMILES IN '+m1+' SET m:matchedMotifProduct_'+str(i)+'_'+str(session_id)+' RETURN (count(m))')
                complete = '"False"'
                if i+1 == len(nodelist):
                    complete = '"True"'
                p2 = session.run('MATCH (s2:matchedSubMotif_'+str(session_id)+')-[r2:motifProductClique]->(p:motifProducts_'+str(i+1)+')<-[r1:motifProductClique]-(s1:matchedMotifProduct_'+str(i)+'_'+str(session_id)+') WHERE p.Complete = '+complete+' AND r1.clique = r2.clique AND ID(r1) <> ID(r2) RETURN ([p.SMILES, r1.clique, s1.SMILES, s2.SMILES, labels(s1), labels(s2)])').value()
                P2 = {}
                for product in p2:
                    p, c, s1, s2, l1, l2 = product
                    n2 = [int(l.split('_')[-1][1:]) for l in l2 if l[:len('matchedSubMotif_'+str(session_id))] == 'matchedSubMotif_'+str(session_id) and len(l) > len('matchedSubMotif_'+str(session_id))]
                    for clique in P1[s1]:
                        n1 = [str(n[0]) for n in clique]
                        for N in n2:
                            if str(N) not in n1:
                                if p not in P2:
                                    P2[p] = []
                                P2[p].extend([tuple(sorted(list(clique)+[(N, s2)], key=lambda x: x[0]))])
                P1 = P2

            for product in P1:
                motifs[product] = list(set(P1[product]))
    
            motifs = {m: list(set(motifs[m])) for m in motifs}
            
            command = []
            for node in nodelist:
                command.extend(['MATCH (n:matchedSubMotif_'+str(session_id)+'_n'+str(node)+') REMOVE n:matchedSubMotif_'+str(session_id)+'_n'+str(node)+' RETURN (count(n))'])
                command.extend(['UNION'])
            for i in range(2, N_ss):
                command.extend(['MATCH (n:matchedMotifProduct_'+str(i)+'_'+str(session_id)+') REMOVE n:matchedMotifProduct_'+str(i)+'_'+str(session_id)+' RETURN (count(n))'])
                command.extend(['UNION'])
            command.extend(['MATCH (n:matchedSubMotif_'+str(session_id)+') REMOVE n:matchedSubMotif_'+str(session_id)+' RETURN (count(n))'])
            session.run('\n'.join(command))

            return [motifs, N_ss]




    def Build(self, csdict, CH2s=None, adjlist=None):
        submotif_matches = {}
        for node in csdict:
            C = list(csdict[node].keys())[0]
            H = csdict[node][C]
            submotif_matches[node] = self.MatchingFunctions.submotifmatcher(C, H)[2]

        no_matches = [] # if any nodes do not match, we return failure
        for node in submotif_matches:
            if len(submotif_matches[node]) == 0:
                no_matches.extend([node])
        if len(no_matches) > 0: # return so user can increase RMSD cutoff
            return 'Nodes',sorted(no_matches,key=lambda x:int(x)),'had no submotif matches...'
        if CH2s: # if user inputted list of CH2 nodes, we can remove CH1 and CH3 matches for that node
            for node in []+CH2s: # in case CH2s outside of the spin-system were inputted, we delete those nodes from the list
                if int(node) not in submotif_matches: # to prevent errors
                    CH2s.remove(node) # means user submitted too much CH2 info since we already ensured all cs nodes matched
            CH2submotifs = self.catalogue_CH2submotifs() # Only keep CH2s
            for node in CH2s:
                node_submotifmatches = submotif_matches[int(node)].copy()
                for submotif in node_submotifmatches:
                    if submotif not in CH2submotifs: 
                        del submotif_matches[int(node)][submotif] # remove CH1 and CH3 for this CH2 node
       
        if len(submotif_matches) != len(csdict):
            return 'Nodes',sorted(no_matches,key=lambda x:int(x)),'had no submotif matches...' # Again, checking to make sure all nodes have matches... CH2 stage could have caused a no match node

        nodes = list(csdict.keys())
        submotifs = []
        for n in submotif_matches:
            submotifs.extend(list(submotif_matches[n].keys()))
        submotifs = set(submotifs)
        valid_pairs = self.NMRT.getvalidpairs()

        neomotifs, N = self.motifBuilderQuery(submotif_matches)

        if N == len(nodes):
            motifs = {}
            for motif in neomotifs:
                if adjlist is not None and len(nodes) > 2:
                    nearest_neighbor = self.nearestneighbor_check(motif, adjlist, neomotifs[motif])
                    if nearest_neighbor is None:
                        continue
                else:
                    nearest_neighbor = neomotifs[motif]
                motifs[motif] = []
                for d in nearest_neighbor:
                    motif_score = sum([submotif_matches[i[0]][i[1]][0]*submotif_matches[i[0]][i[1]][1] for i in d])/float(sum(submotif_matches[i[0]][i[1]][1] for i in d))
                    motifs[motif].extend([(motif, d, round(motif_score,5))])
                motifs[motif] = sorted(motifs[motif], key=lambda x: x[2])
            motifs = sorted([motifs[motif][0] for motif in motifs], key=lambda x: x[2])
            return motifs

        # if need to do some actual building...

        if N != len(nodes):
            # package neomotifs into readable combinations format
            products = {}
            for motif in neomotifs:
                for c in neomotifs[motif]:
                    for combo in itertools.permutations([C[1] for C in c]):
                        if combo not in products:
                            products[combo] = []
                        products[combo].extend([motif])
            for combo in products:
                products[combo] = list(set(products[combo]))
            for combo in itertools.permutations(nodes, N):
                for i in range(N, len(nodes)):
                    choices = [list(submotif_matches[n].keys()) for n in nodes[:i+1]]
                    for c in itertools.product(*choices):
                        if c[:i] in products:
                            pairs = set([tuple(sorted([n, c[i]],key=lambda x:(len(x),x))) for n in c[:i]])
                            query = False
                            for p in pairs:
                                if p in valid_pairs:
                                    query = True
                                    break
                            if query:
                                for pr in products[c[:i]]:
                                    if c[:i+1] not in products:
                                        products[c[:i+1]] = self.adjacency_query(pr, c[i]) 
            motifs = []
            for c in products:
                if len(c) == len(nodes):
                    motifs.extend(products[c])
            return list(set(motifs))

    def newidea(self, spinsys):
        motifproductlist = self.NMRT.assignmotifproducts(N=len(spinsys))   
        Cavg = list(self.MatchingFunctions.csaverages['submotif'][2].keys())[0]
        Havg = self.MatchingFunctions.csaverages['submotif'][2][Cavg]
        motifs = {}
        for product in motifproductlist:
            if motifproductlist[product]['Complete'] == 'False':
                continue
            csdict = copy.deepcopy(motifproductlist[product])
            for node in csdict['replicates']:
                C, c = list(csdict['CS'][node].keys())[0]
                H, h = csdict['CS'][node][(C, c)]
                if csdict['replicates'][node] < self.replicates:
                    csdict['CS'][node] = {(C, Cavg): (H, Havg)}
                if c == 0 and h > 0:
                    csdict['CS'][node] = {(C, Cavg): (H, h)}
                if c > 0 and h == 0:
                    csdict['CS'][node] = {(C, c): (H, Havg)}
                if c == 0 and h == 0:
                    csdict['CS'][node] = {(C, Cavg): (H, Havg)}
            zscore, assignments = self.MatchingFunctions.bipartitescoring(spinsys, csdict, 'zscore')
            if zscore <= self.matching_cutoff:
                motif = product.split('^')[0]
                if motif not in motifs:
                    motifs[motif] = []
                result = [zscore]+sorted([(n, csdict['SMILES'][assignments[n]]) for n in assignments], key=lambda x: x[0])
                motifs[motif].extend([result])
        for motif in motifs:
            motifs[motif] = sorted(motifs[motif], key=lambda x: x[0])[0]
        return motifs



# example spin-systems
#arginine_ss = {8: {43.2: 3.2364}, 4: {57.003: 3.7638}}#, 7: {26.576: 1.6795}, 6: {30.281: 1.909}}

#pm = MotifBuilder()

#pm.Build(csdict=arginine_ss, CH2s=None, adjlist=None)
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
                N = len(matchingDict)
                search_for_completemotifs = True
            else:
                N = N_max
                search_for_completemotifs = False
            #N=2
            nodelist = list(matchingDict.keys())
            # set up session
            command = ['CREATE (Q: MBQ {Session_ID: "'+str(session_id)+'"})'] # create our session
            # add in cs nodes and matching data
            for node in nodelist:
                command.extend(['CREATE (n'+str(node)+':MBQNodes {Node:'+str(node)+', Session_ID:"'+str(session_id)+'"})'])
                command.extend(['CREATE (Q)-[:mbqNode]->(n'+str(node)+')'])
            command.extend(['WITH Q'])
            for node in nodelist:
                command.extend(['MATCH (n:MBQNodes) WHERE n.Node = '+str(node)+' AND n.Session_ID = "'+str(session_id)+'"'])
                command.extend(['MATCH (s:submotif_2) WHERE s.SMILES IN ['+', '.join(['"'+s+'"' for s in matchingDict[node].keys()])+']'])
                command.extend(['CREATE (n)-[:matchedSubMotif]->(s)'])
                command.extend(['SET s:matchedSubMotif_'+str(session_id)+'_n'+str(node)])
                command.extend(['RETURN (n)'])
                command.extend(['UNION'])
            session.run('\n'.join(command[:-1]))
            nodelist = list(matchingDict.keys())
            
            # search for motif products
            motifs = {}          
            if len(nodelist) == N:
                command = ['MATCH (p:motifProducts_'+str(N)+') WHERE p.Complete = "True"']
            else:
                command = ['MATCH (p:motifProducts_'+str(N)+') WHERE p.Complete = "False"']
            for node in nodelist:
                command.extend(['AND (p)<-[:motifProductClique]-(:matchedSubMotif_'+str(session_id)+'_n'+str(node)+')'])
            
            command.extend(['MATCH (p)<-[r:motifProductClique]-(s:submotif_2)']) #
            command.extend(['RETURN ([p.SMILES, r.clique, ID(r), s.SMILES, labels(s)])'])
            motif_hits = session.run('\n'.join(command)).value()

            cliques = {M:{} for M in set([m[0] for m in motif_hits])}
            for m in motif_hits:
                if len(m[4]) > 1:
                    if m[1] not in cliques[m[0]]:
                        cliques[m[0]][m[1]] = {}
                    for s in m[4]:
                        if s[:len('matchedSubMotif_'+str(session_id))] == 'matchedSubMotif_'+str(session_id):
                            node = int(s[len('matchedSubMotif_'+str(session_id))+2:])
                            if node not in cliques[m[0]][m[1]]:
                                cliques[m[0]][m[1]][node] = []
                            cliques[m[0]][m[1]][node].extend([(m[3], m[2])])
            # cliques: motif: clique_number: node: [(submotif smiles, relationship ID)]
            motifs = {}
            for motif in cliques:
                for clique in cliques[motif]:
                    if len(cliques[motif][clique]) != N:
                        continue
                    nodes = list(cliques[motif][clique].keys())
                    submotifs = [cliques[motif][clique][node] for node in nodes]
                    for combo in itertools.product(*submotifs):
                        rel_ids = [c[1] for c in combo]
                        if len(set(rel_ids)) != N:
                            continue
                        if motif not in motifs:
                            motifs[motif] = []
                        motifs[motif].extend([tuple(sorted([(nodes[i], combo[i][0]) for i in range(len(nodes))], key=lambda x: x[0]))])  
            #print(session.run('MATCH (n) RETURN distinct labels(n)').value())
            command = []
            for node in nodelist:
                command.extend(['MATCH (n:matchedSubMotif_'+str(session_id)+'_n'+str(node)+') REMOVE n:matchedSubMotif_'+str(session_id)+'_n'+str(node)+' RETURN (n)'])
                command.extend(['UNION'])
            session.run('\n'.join(command[:-1]))

            return [motifs, N]




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
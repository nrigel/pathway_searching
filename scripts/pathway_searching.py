from neo4j_connect import NMRTServer, ReactomeServer
from matching_functions import MatchingFunctions

class PathwayMapper(object):

    def __init__(self, matching_cutoff=5, replicates=3):
        # store our initial variables
        self.matching_cutoff, self.replicates = matching_cutoff, replicates
        self.r = (1, 5)
        
        self.NMRT, self.Reactome = NMRTServer(), ReactomeServer()

        self.MatchingFunctions = MatchingFunctions(self.NMRT, self.matching_cutoff, self.replicates)
            
    def reactomequery(self, submotifmatches, motifmatches, metabolitematches):
        pathway_matches = {'submotifs': {}, 'motifs': {}, 'metabolites': {}}

        self.Reactome.generategraphs() # generate a few networkx subgraphs
        
        # return pathway match list for each node/motif
        for metabolite in metabolitematches:
            pathway_matches['metabolites'][metabolite] = self.Reactome.pathwaymatch(metabolitematches[metabolite]['SMILES'], structuretype='metabolites', shell=None)

        for node in submotifmatches:
            pathway_matches['submotifs'][node] = {}
            for shell in submotifmatches[node]:
                pathway_matches['submotifs'][node][shell] = {}
                for submotif in submotifmatches[node][shell]:
                    pathway_matches['submotifs'][node][shell][submotif] = self.Reactome.pathwaymatch(submotif, structuretype='submotifs', shell=shell)

        for shell in motifmatches:
            pathway_matches['motifs'][shell] = {}
            for motif in motifmatches[shell]:
                # Reactome does not include spin-systems in the SMILES
                pathway_matches['motifs'][shell][motif] = self.Reactome.pathwaymatch(motif.split('_')[0], structuretype='motifs', shell=shell)
                # motif: {species: pathway list}
        return pathway_matches

    def spinsystem(self, ss, motiflist={}):
        metabolite_matches = self.MatchingFunctions.metabolitematcher(ss)
        submotif_matches = {}
        for node in ss:
            for C in ss[node]:
                submotif_matches[node] = self.MatchingFunctions.submotifmatcher(C, ss[node][C])
        motif_matches = self.MatchingFunctions.motifmatcher(ss)    
        for motif in motiflist:
            motif_matches[2][motif] = motiflist[motif]
        pathway_matches = self.reactomequery(submotif_matches, motif_matches, metabolite_matches)
        return pathway_matches


        # next steps
            # integrate motif builder/write smaller script than we have and call it from here
            # check that pathway matching to the motiflist inputted above actually works
            # update reactome database?
                # add metabolites that are not in our databases
                # review what we did last time
                # make new one with all possible compounds
                # then make sub-motif and motif portions that we can search for
                # should be able to build database without chemical shifts
            # write reactome function 
            # write short web-server to demonstrate what we can do
                # include example amino acids with maybe even spin-system graphs
                # drop down menu: select example compound, match for sub-motifs, motifs and search those against Reactome

            # instead of clique matching why not just create a product node for each clique with attached sub-motif chemical shifts and then do scoring against those??

    
# example spin-systems
#arginine_ss = {8: {43.2: 3.2364}, 4: {57.003: 3.7638}, 7: {26.576: 1.6795}, 6: {30.281: 1.909}}

#pm = PathwayMapper()
#pm.spinsystem(arginine_ss)

#NMRT.close()
#Reactome.close()
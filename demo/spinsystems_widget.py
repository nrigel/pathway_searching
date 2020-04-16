import sys
sys.path.append('../scripts')
from pathway_searching import PathwayMapper
from motif_builder import MotifBuilder

from PySide2.QtWidgets import QWidget

class SpinSystems(QWidget):

    def __init__(self):
        QWidget.__init__(self) 

        self.pathwaymapper = PathwayMapper()
        self.spinsystems = self.retrieveSpinSystemList()
        self.MotifBuilder = MotifBuilder()
        self.built_motifs = {}

    def retrieveSpinSystemList(self):
        spinsystems = {}
        metabolites = self.pathwaymapper.NMRT.metabolitegraph()
        for metabolite in metabolites['spinsystems']:
            for spinsys in metabolites['spinsystems'][metabolite]:
                    csdict = {}
                    for node in spinsys.split():
                        neo = metabolites['networkx'].nodes[metabolite+'_'+node]['neo4j']
                        if 'CS' in neo.keys():
                            csdict[int(node)] = {float(neo['CS'].split()[0]): float(neo['CS'].split()[1])}
                    if len(csdict) == len(spinsys.split()):
                        spinsystems[metabolite+'_'+spinsys] = csdict
        return spinsystems

    def motifbuilder(self, spinsystem):
        motifs = self.MotifBuilder.Build(csdict=self.spinsystems[spinsystem], CH2s=None, adjlist=None)
        # matches[i][motif] = {'zscore': zscore, 'assignments': assignments, 'CS': csdict}
        self.built_motifs = {}
        for motif in motifs:
            if motifs[motif][-1] <= self.MotifBuilder.matching_cutoff:
                self.build_motifs[motif] = motifs[motif]

    def displayspinsystem(self, spinsystem):
        # show some stuff
        cstext = ['Node     C   H']
        for node in sorted(list(self.spinsystems[spinsystem].keys()), key=lambda x:x):
            C = list(self.spinsystems[spinsystem][node].keys())[0]
            H = self.spinsystems[spinsystem][node][C]
            cstext.extend([str(node)+'  '+str(C)+'  '+str(H)])
        self.main_widget.csbox.setText('\n'.join(cstext))

    def pathwaysearch(self, spinsystem):
        # search for matching pathways
        self.pathwaymatches = self.pathwaymapper.spinsystem(self.spinsystems[spinsystem], self.built_motifs)
        self.showpathwaylist()

    def showpathwaylist(self):
        opts = [k for k in self.main_widget.structure_menu.keys() if k != 'layout' and self.main_widget.structure_menu[k].isChecked()]
        # show pathway list
        pathwaylist = {}
        # metabolites -> metabolite -> species -> pathway list
        if 'metabolites' in opts:
            for metabolite in self.pathwaymatches['metabolites']:
                for species in self.pathwaymatches['metabolites'][metabolite]:
                    if species not in pathwaylist:
                        pathwaylist[species] = set()
                    for pathway in self.pathwaymatches['metabolites'][metabolite][species]:
                        pathwaylist[species].add(pathway)
        # motifs -> shell -> species -> pathway list
        for shell in self.pathwaymatches['motifs']:
            if 'motif_'+str(shell) not in opts:
                continue
            for motif in self.pathwaymatches['motifs'][shell]:
                for species in self.pathwaymatches['motifs'][shell][motif]:
                    if species not in pathwaylist:
                        pathwaylist[species] = set()
                    for pathway in self.pathwaymatches['motifs'][shell][motif][species]:
                        pathwaylist[species].add(pathway)
        # submotifs -> node -> shell -> species -> pathway list
        for node in self.pathwaymatches['submotifs']:
            for shell in self.pathwaymatches['submotifs'][node]:
                if 'submotif_'+str(shell) not in opts:
                    continue
                for submotif in self.pathwaymatches['submotifs'][node][shell]:
                    for species in self.pathwaymatches['submotifs'][node][shell][submotif]:
                        if species not in pathwaylist:
                            pathwaylist[species] = set()
                        for pathway in self.pathwaymatches['submotifs'][node][shell][submotif][species]:
                            pathwaylist[species].add(pathway)
        
        self.main_widget.pathwaylist.clear()
        if self.main_widget.species_menu.currentText() in pathwaylist:
            self.main_widget.pathwaylist.addItems(sorted(list(pathwaylist[self.main_widget.species_menu.currentText()]), key=lambda x: x))
        self.main_widget.pathwaycount.setText('Count: '+str(self.main_widget.pathwaylist.count()))

        G = self.pathwaymapper.Reactome.graphs['pathways_to_processes']['networkx']
        nodes = {}
        for node in G.nodes:
            if G.nodes[node]['name'] not in nodes:
                nodes[G.nodes[node]['name']] = []
            nodes[G.nodes[node]['name']].extend([node])
        processeslist = set()
        for pathway in pathwaylist[self.main_widget.species_menu.currentText()]:
            if pathway in nodes:
                for p in nodes[pathway]:
                    for edge in G.edges(p):
                        for e in edge:
                            if e != p:
                                processeslist.add(G.nodes[e]['name'])
        self.main_widget.processlist.addItems(sorted(list(processeslist), key=lambda x: x))
        self.main_widget.processcount.setText('Count: '+str(self.main_widget.processlist.count()))
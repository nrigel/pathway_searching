import openbabel, pybel, time
import networkx as nx
import networkx.algorithms.isomorphism as iso
#from rdkit import Chem
#from rdkit.Chem.Draw import rdMolDraw2D
#https://baoilleach.blogspot.com/2010/11/automorphisms-isomorphisms-symmetry.html

class Compound(object): # my class for molecules, metabolites, motifs, submotifs, etc.
    def __init__(self, structure):
        self.mapping = None
        if type(structure) == str:
            self.smiles = structure # SMILES is our input, either canonical or not
            self.obmol = self.smiles_to_obmol() # convert to an OBMol object
            self.G = self.obmol_to_graph() # translate OBMol to nxGraph for easier handling
            #self.pybelmol = self.obmol_to_pybelmol() # keep pybelMol object for requests
        if type(structure) == openbabel.OBMol:
            self.obmol = structure # OBMol object is the input
            self.G = self.obmol_to_graph()
            self.smiles, self.canonicalOrder = self.obmol_to_smiles()
            #self.pybelmol = self.obmol_to_pybelmol()
        if type(structure) == pybel.Molecule:
            self.pybelmol = structure # PybelMol object is the input
            self.obmol = self.pybelmol_to_obmol()
            self.G = self.obmol_to_graph()
            self.smiles, self.canonicalOrder = self.obmol_to_smiles()
        if type(structure) == nx.classes.graph.Graph:
            self.G = structure.copy() # nxGraph is the input, we will not renumber the atoms in the graph, will update self.mapping
            self.obmol = self.G_to_obmol()
            self.smiles, self.canonicalOrder = self.obmol_to_smiles()
            #self.pybelmol = self.obmol_to_pybelmol()
    
    def update(self): # changes to G require update to have proper MOL
        self.obmol = self.G_to_obmol() # remake OBMol object with current nxGraph
        self.smiles, self.canonicalOrder = self.obmol_to_smiles() # rewrite SMILES
        #self.pybelmol = self.obmol_to_pybelmol() # remake pybelmol

    def canonicalize(self):
        start = time.time()
        self.check_charges()
        mapping = {}
        if self.mapping:
            aliases = {self.mapping[n]: n for n in self.mapping}
        else:
            aliases = {n: n for n in list(self.G.nodes)}
        smiles, order = self.obmol_to_smiles()
        order = order.split() # node indices in order of new smiles string
        for node in order:
            mapping[aliases[int(node)]] = order.index(node)+1
        S = Compound(smiles)
        self.G = S.G
        self.update()
        return mapping
    
    def obmol_to_pybelmol(self):
        return pybel.Molecule(self.obmol) # convert OBMol object to pybelmol object
    
    def pybelmol_to_obmol(self):
        obConversion = openbabel.OBConversion() # my strategy is to make the MOL table and then convert that into OBMol
        obConversion.SetInAndOutFormats('mol', 'mol')
        ob_mol = openbabel.OBMol()
        obConversion.ReadString(ob_mol, self.pybelmol.write('mol'))
        return ob_mol
    
    def smiles_to_obmol(self):
        obConversion = openbabel.OBConversion() # input the smiles into OBMol object, check indices maybe
        obConversion.SetInAndOutFormats('smi', 'smi')
        ob_mol = openbabel.OBMol()
        obConversion.ReadString(ob_mol, self.smiles)
        ob_mol.DeleteHydrogens()
        return ob_mol
    
    def obmol_to_smiles(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats('mol', 'smi')
        obConversion.AddOption("c", obConversion.OUTOPTIONS) # set to canonical
        obConversion.AddOption("i", obConversion.OUTOPTIONS) # ignore chirality
        smiles = obConversion.WriteString(self.obmol)[:-2]
        atomOrder = self.obmol.GetData("SMILES Atom Order").GetValue()
        return smiles, atomOrder
    
    def simplify_obmol(self,negcharges=False):
        self.check_charges(negcharges) # remove pos charge artifacts and all negative charges
        self.update() # reset OBMol, SMILES, and pybelMol
        return self.smiles # return the updated canonical smiles

    def check_charges(self,negcharges=False): # negcharges True means that they are allowed
        atoms = list(openbabel.OBMolAtomIter(self.obmol))
        for atom in self.G.nodes:
            a = atom-1 # G indices are OBMol indices+1
            if self.mapping: # means OBMol was loaded from G at some point
                a = self.mapping[atom]-1 # G indices are OBMol indices+1
            if self.G.nodes[atom]['Charge'] > 0: # iterate thru charged atoms
                if atoms[a].ExplicitHydrogenCount(): # check if OBMol bonded a hydrogen to the atom to compensate for false charge
                    self.G.nodes[atom]['Charge'] = 0 # reset charge in self.G
                    self.obmol.BeginModify() # need this to make changes to existing OBMol object
                    atoms[a].SetFormalCharge(0) # delete charge (yes this will cause problems for >+1, but best method I could come up with)
                    self.obmol.EndModify() # need this to save changes to existing OBMol object
                else: # set charge to 0, check valences, re-adjust
                    if atoms[a].ImplicitHydrogenCount() > 0:
                        self.G.nodes[atom]['Charge'] = 0
                        self.obmol.BeginModify()
                        atoms[a].SetFormalCharge(0)
                        self.obmol.EndModify()
            if self.G.nodes[atom]['Charge'] < 0 and negcharges is False: # means we will remove this neg charge
                self.G.nodes[atom]['Charge'] = 0 # set to 0
        self.remove_hydrogens() # remove hydrogens that aren't needed anymore
        self.update() # reset OBMol, SMILES, and pybelMol
        return 
    
    def check_valences(self,node=None,add_charges=True):
        if node:
            nodes = [node]
        else:
            nodes = list(self.G.nodes)
        self.remove_hydrogens() # remove hydrogens that aren't needed anymore
        incorrect,valences = [],{6:4,7:3,8:2}
        for atom in nodes:
            bondorder = self.get_hvybondorder(atom)
            if self.G.nodes[atom]['AtomicNum'] in valences:
                if bondorder > valences[self.G.nodes[atom]['AtomicNum']] and self.G.nodes[atom]['Charge'] == 0:
                    if add_charges and self.G.nodes[atom]['AtomicNum'] != 6:
                        charge = bondorder - valences[self.G.nodes[atom]['AtomicNum']]
                        if charge == 1: # add pos charge
                            self.G.nodes[atom]['Charge'] = charge
                        elif charge == 0:
                            continue
                        else:
                            incorrect.extend([str(atom)])
                    else:
                        incorrect.extend([str(atom)])
        self.update() # reset OBMol, SMILES, and pybelMol
        return incorrect
    
    def get_hvybondorder(self,node):
        bondorder = 0
        for n in self.D()[node]:
            if self.G.nodes[n]['AtomicNum'] == 1: # skip hydrogens
                continue 
            if (node,n) in self.G.edges: # should be since it is the adjacency matrix but here just in case
                bondorder += self.G.edges[(node,n)]['BondOrder'] # add bond to running bond order
        return bondorder
    
    def reset_aromaticity(self): # needed to fix BO values of generated motifs and sub-motifs   
        if self.mapping:
            rmapping = {self.mapping[A]: A for A in self.mapping}
        for bond in openbabel.OBMolBondIter(self.obmol):
            if bond.IsAromatic():
                BO = 4
            else:
                BO = bond.GetBondOrder()
            a1,a2 = bond.GetBeginAtom().GetIndex()+1,bond.GetEndAtom().GetIndex()+1
            if self.mapping:
                if a1 in rmapping and a2 in rmapping:
                    a1,a2 = rmapping[a1],rmapping[a2]
                else:
                    a1,a2 = 0,0
            if (a1,a2) in self.G.edges:
                self.G.edges[(a1,a2)]['BO'] = BO
        return self.smiles
    
    def descriptors(self,option=None):
        self.pybelmol = self.obmol_to_pybelmol()
        if option == 'logP':
            return self.pybelmol.calcdesc()['logP']
    
    def MOL(self): # write MOL output for quick 2D visualization
        mymol = pybel.Molecule(self.obmol)
        mymol.make3D() # generate 3D coordinates to avoid error message printed in terminal
        mymol.removeh() # remove hydrogens
        mymol.draw(show=False, filename=None, update=True, usecoords=False) # generate 2D coordinates
        return mymol.write('mol')
    
    def Draw(self,spinsystem=None,firstshell=[],secondshell=[],labels={},molSize=(450,150)):
        # first need to check for explicit hydrogens because RDKit ignores them
        if 'H' in self.smiles:
            # need to remove it
            S = Compound(self.G.copy())
            S.remove_hydrogens()
            aliases = {}
            for node in S.G.nodes:
                aliases[node] = len(aliases)+1
            MOL = S.MOL()
        else:
            MOL = self.MOL()
            aliases = {node:node for node in self.G}
        if len(labels) != 0:
            labels = {aliases[n]-1:labels[n] for n in labels}
        if spinsystem is None: # generate atom details if needed
            spinsystem = []
            for ss in self.spinsystems():
                spinsystem.extend(ss.split())
            
            firstshell,secondshell = [],[]
            for atom in spinsystem:
                for node in self.D()[int(atom)]:
                    if str(node) not in spinsystem:
                        firstshell.extend([str(node)])
            for atom in firstshell:
                for node in self.D()[int(atom)]:
                    if str(node) not in spinsystem and str(node) not in firstshell:
                        secondshell.extend([str(node)])
            spinsystem = list(set([aliases[int(a)] for a in spinsystem]))
            firstshell = list(set([aliases[int(a)] for a in firstshell]))
            secondshell = list(set([aliases[int(a)] for a in secondshell]))
            if len(labels) == 0:
                labels = {}
                for atom in spinsystem:
                    labels[atom-1] = 'C'+str(atom)
        else:
            spinsystem = [aliases[int(n)] for n in spinsystem]
            firstshell = [aliases[int(n)] for n in firstshell]
            secondshell = [aliases[int(n)] for n in secondshell]
        # colors are RGB but normalized (/255)
        color_scheme = {'spinsystem':(0.91,0.4,0.4),'firstshell':(0.4,0.55,0.91),'secondshell':(0.4,0.86,0.66)}
        colors,highlight = {},[]
        for atom in spinsystem:
            colors[atom-1] = color_scheme['spinsystem']
            highlight.extend([atom-1])
        for atom in firstshell:
            colors[atom-1] = color_scheme['firstshell']
            highlight.extend([atom-1])
        for atom in secondshell:
            colors[atom-1] = color_scheme['secondshell']
            highlight.extend([atom-1])

        # http://rdkit.blogspot.com/2015/02/new-drawing-code.html
        
        mol = Chem.MolFromMolBlock(MOL)
        mc = Chem.Mol(mol.ToBinary())
        Chem.rdDepictor.Compute2DCoords(mc)
        remove_aromatics = False # remove aromatic atoms -- doesnt work yet
        if remove_aromatics:
            # load coordinates from RDKit into OpenBabel to remove aromaticity
            obConversion = openbabel.OBConversion() # input the smiles into OBMol object, check indices maybe
            obConversion.SetInAndOutFormats('mol', 'mol')
            ob_mol = openbabel.OBMol()
            obConversion.ReadString(ob_mol, Chem.MolToMolBlock(mc))
            # remove aromatic properties
            ob_mol.UnsetAromaticPerceived()
            for atom in openbabel.OBMolAtomIter(ob_mol):    
                if atom.IsAromatic():
                    atom.UnsetAromatic()
            for bond in openbabel.OBMolBondIter(ob_mol):
                if bond.IsAromatic():
                    bond.UnsetAromatic()
            # use pybel to generate MOL block
            mymol = pybel.Molecule(ob_mol)
            unset_aromaticmol = mymol.write('mol')
            print(mymol.write('smi'))
            # load new mol block back into rdkit
            mol = Chem.MolFromMolBlock(unset_aromaticmol)
            mc = Chem.Mol(mol.ToBinary())

        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
        opts = drawer.drawOptions()
        for atom in labels:
            opts.atomLabels[atom] = labels[atom]
        drawer.DrawMolecule(mc,highlightAtoms=highlight,highlightAtomColors=colors,highlightBonds=[])
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace('svg:','')
        # need to correct the y-axis cutoff in the svg with <g transform="translate(0,15)">
        svg = svg.split('\n') # break into list for editing purposes
        # correct the height
        whline = svg[6].split() # call width, height line
        hchange = 30
        newheight = int(whline[1][8:whline[1].find('p')])+hchange # add 30 points to current height
        whline[1] = whline[1][:8]+str(newheight)+whline[1][whline[1].find('p'):] # change height portion to new height
        svg[6] = ' '.join(whline) # change width-height line to use new height
        # correct rectangle height, location, and opacity
        opacity = 1 # 0 makes background transparent, 1 makes background white/whatever you set rectangle to
        rect = svg[7].split() # <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='450' height='180' x='0' y='-15'> </rect>
        opacline = rect[1].split(':')
        opacline[1] = str(opacity)+opacline[1][opacline[1].find(';'):]
        rect[1] = ':'.join(opacline)
        rect[3] = rect[3][:7]+"'"+str(newheight)+"'"
        rect[5] = rect[5][:2]+"'-"+str(int(hchange/2.0))+"'>"
        svg[7] = ' '.join(rect)
        # now we need to add in an object delimiter so we can shift this whole thing down
        svg.insert(7,"<g transform='translate(0,15)'>")
        # now we need to add indents to every line after 7
        for l in range(8,len(svg)-2):
            svg[l] = '  '+svg[l]
        # now we add in the line to end the object
        svg.insert(-2,"</g>")
        return '\n'.join(svg) # rejoin svg

    def obmol_to_graph(self): # convert obmol object to nx graph
        G = nx.Graph()
        atoms = list(openbabel.OBMolAtomIter(self.obmol))
        for a in range(len(atoms)):
            G.add_node(a+1,AtomicNum=atoms[a].GetAtomicNum(),Charge=atoms[a].GetFormalCharge())
        for bond in openbabel.OBMolBondIter(self.obmol):
            if bond.IsAromatic():
                BO = 4
            else:
                BO = bond.GetBondOrder()
            G.add_edge(bond.GetBeginAtom().GetIndex()+1,bond.GetEndAtom().GetIndex()+1,BO=BO,BondOrder=bond.GetBondOrder())
        return G
    
    def G_to_obmol(self): # graph to pybel mol
        node_data = self.G.nodes(data=True)
        node_data = sorted(node_data,key=lambda x:x[0])
        mymol = openbabel.OBMol()
        atoms = {}
        for n in range(len(node_data)):
            node = node_data[n]
            a = mymol.NewAtom()
            a.SetAtomicNum(node[1]['AtomicNum'])
            a.SetFormalCharge(node[1]['Charge'])
            atoms[node[0]] = n+1
        for edge in self.G.edges(data=True):
            mymol.AddBond(atoms[edge[0]],atoms[edge[1]],edge[2]['BondOrder'])
            if edge[2]['BO'] == 4:
                mymol.GetBond(atoms[edge[0]],atoms[edge[1]]).SetAromatic()
                mymol.GetAtom(atoms[edge[0]]).SetAromatic()
                mymol.GetAtom(atoms[edge[1]]).SetAromatic()
        self.mapping = atoms
        return mymol
    
    # some graph functions
    
    def D(self): # generate adjacency lists for G without needing the full command
        return nx.to_dict_of_lists(self.G)
    
    def get_adjlist(self, nodes):
        nodes = [str(n) for n in nodes]
        adjlist = []
        D = self.D()
        for node in nodes:
            a = []
            for n in D[int(node)]:
                if str(n) in nodes:
                    a.extend([int(n)])
            adjlist.extend([tuple([int(node), tuple(sorted(a, key=lambda x: x))])])
        return sorted(adjlist, key=lambda x: x[0])

    def remove_hydrogens(self): # remove hydrogens if they exist, make isomorphisms fail if present
        H_count = 0
        nodes = list(self.G.nodes)
        for node in nodes:
            if self.G.nodes[node]['AtomicNum'] == 1:
                self.G.remove_node(node)
                H_count += 1
        return str(H_count)+' hydrogens removed.'
    
    def cyclic_nodes(self): # generate list of nodes that are in rings
        cycles = list(nx.cycle_basis(self.G))
        cyclic_nodes = []
        for cycle in cycles:
            str_cycle = [str(node) for node in cycle]
            cyclic_nodes.extend(str_cycle)
        cyclic_nodes = sorted(list(set(cyclic_nodes)),key=lambda x:int(x))
        return cyclic_nodes 
    
    def CHs(self,NODE=None): # generate list of nodes that are CH2s
        if NODE:
            bondorder = self.get_hvybondorder(node=NODE)
            if self.get_hvybondorder(node=NODE) == 3:
                return True
            else:
                return
        chs = []
        for node in self.G.nodes:
            if self.G.nodes[node]['AtomicNum'] != 6:
                continue
            if self.get_hvybondorder(node=node) == 3:
                chs.extend([str(node)])
        return sorted(chs,key=lambda x:int(x))
    
    def CH2s(self,NODE=None): # generate list of nodes that are CH2s
        if NODE:
            bondorder = self.get_hvybondorder(node=NODE)
            if self.get_hvybondorder(node=NODE) == 2:
                return True
            else:
                return
        ch2s = []
        for node in self.G.nodes:
            if self.G.nodes[node]['AtomicNum'] != 6:
                continue
            if self.get_hvybondorder(node=node) == 2:
                ch2s.extend([str(node)])
        return sorted(ch2s,key=lambda x:int(x))

    def CH3s(self,NODE=None): # generate list of nodes that are CH2s
        if NODE:
            bondorder = self.get_hvybondorder(node=NODE)
            if self.get_hvybondorder(node=NODE) == 1:
                return True
            else:
                return
        ch3s = []
        for node in self.G.nodes:
            if self.G.nodes[node]['AtomicNum'] != 6:
                continue
            if self.get_hvybondorder(node=node) == 1:
                ch3s.extend([str(node)])
        return sorted(ch3s,key=lambda x:int(x))
    
    def connected_component_subgraphs(self, G=None):
        if G is None:
            G = self.G
        subgraphs = []
        for c in nx.connected_components(G):
            subgraphs.extend([G.subgraph(c)])
        return subgraphs

    def spinsystems(self,NODE=None): # generate list of spin-systems
        H_count = {}
        if len(self.G.edges) == 0 and self.G.nodes[list(self.G.nodes)[0]]['AtomicNum'] == 6:
            H_count[list(self.G.nodes)[0]] = 4 # means it is a single CH4
        for edge in self.G.edges(data=True):
            if self.G.nodes[edge[0]]['AtomicNum'] == 1 or self.G.nodes[edge[1]]['AtomicNum'] == 1:
                continue
            for node in edge[:2]:
                if self.G.nodes[node]['AtomicNum'] == 6:
                    if node not in H_count:
                        H_count[node] = 4
                    H_count[node] -= edge[2]['BondOrder']
        H = self.G.copy()
        for node in self.G.nodes:
            if self.G.nodes[node]['AtomicNum'] != 6 or H_count[node] == 0:
                H.remove_node(node)
        spinsystems = self.connected_component_subgraphs(H)
        spinsystems = [list(spinsys.nodes) for spinsys in spinsystems]
        for s in range(len(spinsystems)):
            spinsystems[s] = ' '.join([str(n) for n in sorted(spinsystems[s],key=lambda x:x)])
        if NODE is None:
            return spinsystems
        for spinsys in spinsystems:
            if str(NODE) in spinsys.split():
                return spinsys
    
    def motifs(self,shell): # generate motifs for the molecule with desired shell (0,1,2 shells only)
        motif_graphs = {}
        for spinsys in self.spinsystems():
            ss_G = nx.Graph()
            ss = spinsys.split()
            for node in ss:
                if node not in ss_G.nodes:
                    ss_G.add_node(int(node),AtomicNum=self.G.nodes[int(node)]['AtomicNum'],Charge=self.G.nodes[int(node)]['Charge'])
                for n in self.D()[int(node)]:
                    if str(n) in ss or shell > 0:
                        edge_data = self.G.get_edge_data(int(node),n)
                        if n not in ss_G.nodes:
                            ss_G.add_node(n,AtomicNum=self.G.nodes[n]['AtomicNum'],Charge=self.G.nodes[n]['Charge'])
                        ss_G.add_edge(int(node),n,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
                    if shell > 1:
                        for N in self.D()[n]:
                            edge_data = self.G.get_edge_data(n,N)
                            if N not in ss_G.nodes:
                                ss_G.add_node(N,AtomicNum=self.G.nodes[N]['AtomicNum'],Charge=self.G.nodes[N]['Charge'])
                            ss_G.add_edge(n,N,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
            motif_graphs[spinsys] = ss_G
        return motif_graphs    
    
    def sub_motifs(self,node): # generates sub-motifs 1,2,3,4 shells -- could even go farther if desired
        sub_motif_graphs = {}
        for i in range(1,5):
            sub_motif_graphs[i] = nx.Graph()
            I,nodes = 0,[node]
            sub_motif_graphs[i].add_node(node,AtomicNum=self.G.nodes[node]['AtomicNum'],Charge=self.G.nodes[node]['Charge'])
            while i != I:
                I += 1
                NODES = []
                for n in nodes:
                    for N in self.D()[n]:
                        NODES.extend([N])
                        edge_data = self.G.get_edge_data(n,N)
                        if n not in sub_motif_graphs[i].nodes:
                            sub_motif_graphs[i].add_node(n,AtomicNum=self.G.nodes[n]['AtomicNum'],Charge=self.G.nodes[n]['Charge'])
                        if N not in sub_motif_graphs[i].nodes:
                            sub_motif_graphs[i].add_node(N,AtomicNum=self.G.nodes[N]['AtomicNum'],Charge=self.G.nodes[N]['Charge'])
                        sub_motif_graphs[i].add_edge(n,N,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
                nodes = []+NODES
        return sub_motif_graphs

    def sub_motif_smiles(self,node):
        sub_motif_graphs = self.sub_motifs(node)
        smiles = {}
        for shell in sub_motif_graphs:
            s = Compound(sub_motif_graphs[shell])
            m = s.canonicalize()
            smiles[shell] = s.smiles+'_'+str(m[node])
        return smiles
    
    def submotif_fragmenter(self,node,shell=2): # generate 2nd shell motif fragments for motif building 
        fragments = {node:self.G.copy()}
        for n in self.D()[node]:
            fragments[n] = nx.Graph()
            fragments[n].add_node(node,AtomicNum=self.G.nodes[node]['AtomicNum'],Charge=self.G.nodes[node]['Charge'])
            edge_data = self.G.get_edge_data(node,n)
            fragments[n].add_node(n,AtomicNum=self.G.nodes[n]['AtomicNum'],Charge=self.G.nodes[n]['Charge'])
            fragments[n].add_edge(node,n,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
            
            if shell > 1:
                for N in self.D()[node]:
                    if n == N:
                        continue
                    edge_data = self.G.get_edge_data(node,N)
                    fragments[n].add_node(N,AtomicNum=self.G.nodes[N]['AtomicNum'],Charge=self.G.nodes[N]['Charge'])
                    fragments[n].add_edge(node,N,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
                for N in self.D()[n]:
                    edge_data = self.G.get_edge_data(n,N)
                    if N not in fragments[n].nodes:
                        fragments[n].add_node(N,AtomicNum=self.G.nodes[N]['AtomicNum'],Charge=self.G.nodes[N]['Charge'])
                    fragments[n].add_edge(n,N,BondOrder=edge_data['BondOrder'],BO=edge_data['BO'])
        del fragments[node]
        for n in fragments:
            fragments[n] = Compound(fragments[n]) # ~0.001 seconds
            SG = dict([[v,k] for k,v in fragments[n].canonicalize().items()]) # ~0.001 seconds
            fragments[n].G = nx.relabel_nodes(fragments[n].G,SG) # ~0.0001 seconds
            fragments[n].update() # ~0.0002 seconds
        return fragments

    def get_termini(self,spinsys=None): # generate the end carbons of linear spin-systems for motif enhancing
        if spinsys is None:
            spinsystems = self.spinsystems() # untargeted
        else:
            spinsystems = [spinsys]
        termini = [] # spin system atoms connected to heteroatoms
        for spinsystem in spinsystems:
            for node in self.D():
                if str(node) not in spinsystem.split():
                    continue
                hetatoms = self.get_heteroatoms(spinsystem)
                for n in self.D()[node]:
                    if str(n) in hetatoms: # connected to heteroatom (makes it a terminus in the spin system)
                        termini.extend([str(node)])
                        break
        return termini
    
    def get_heteroatoms(self,spinsys=None): # make list of atoms not in spinsystem
        if spinsys is None:
            spinsystem = []
            for ss in self.spinsystems(): # untargeted
                spinsystem.extend(ss.split())
        else:
            spinsystem = spinsys.split()
        hetatoms = []
        for node in self.G.nodes:
            if str(node) not in spinsystem:
                hetatoms.extend([str(node)])
        return sorted(list(set(hetatoms)),key=lambda x: int(x)) 

    def share_spinsystem_check(self,nodes):
        spinsystems = []
        for spinsystem in self.spinsystems():
            common = [n for n in spinsystem.split() if n in [str(N) for N in nodes]]
            if len(common) > 0:
                spinsystems.extend([(spinsystem)])
        return spinsystems



    

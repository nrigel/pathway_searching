from PySide2.QtCore import Qt
from PySide2.QtSvg import QSvgWidget
from PySide2.QtWidgets import QGridLayout, QHBoxLayout, QVBoxLayout, QTabWidget, QWidget, QComboBox, QLabel, QCheckBox, QListWidget, QPushButton
from spinsystems_widget import SpinSystems

class MainWidget(QWidget):
    
    def __init__(self):
        QWidget.__init__(self)  
        
        #self.grabKeyboard()

    def initialize(self):

        self.spinsystems_widget = SpinSystems()
        self.spinsystems_widget.main_widget = self

        # set up main layout
        self.mainLayout()
        
    def mainLayout(self):
        self.main_layout = QHBoxLayout()
        self.left_layout = QVBoxLayout()
        self.main_layout.addLayout(self.left_layout)
        self.main_layout.setStretchFactor(self.left_layout, 2)
        self.setLayout(self.main_layout)
        
        self.dropdownmenu()
        self.structure_disp = QSvgWidget()
        self.csbox = QLabel()
        self.spinsystems_widget.displayspinsystem(self.spinsystem_menu.currentText())
        self.left_layout.addWidget(self.csbox)

        self.left_layout.addWidget(self.structure_disp)

        self.motif_builder = QPushButton('Build Motifs')
        def buildmotifs():
            self.spinsystems_widget.motifbuilder(self.spinsystem_menu.currentText())
        self.motif_builder.clicked.connect(buildmotifs)
        self.left_layout.addWidget(self.motif_builder)
        
        self.map_pathways = QPushButton('Map Pathways')
        def mappathways():
            self.spinsystems_widget.pathwaysearch(self.spinsystem_menu.currentText())
            self.spinsystems_widget.showpathwaylist()
        self.map_pathways.clicked.connect(mappathways)
        self.left_layout.addWidget(self.map_pathways)

        self.structureoptions()
        self.speciesmenu()

        self.builtmotiflayout = QVBoxLayout()
        self.builtmotiflist = QListWidget()
        self.builtmotiflayout.addWidget(self.builtmotiflist)
        self.builtmotiflayout.setStretchFactor(self.builtmotiflist, 9)
        self.builtmotifcount = QLabel('Count: ')
        self.builtmotiflayout.addWidget(self.builtmotifcount)
        self.builtmotiflayout.setStretchFactor(self.builtmotifcount, 1)
        self.main_layout.addLayout(self.builtmotiflayout)
        self.main_layout.setStretchFactor(self.builtmotiflayout, 4)

        self.pathwaylistlayout = QVBoxLayout()
        self.pathwaylist = QListWidget()
        self.pathwaylistlayout.addWidget(self.pathwaylist)
        self.pathwaylistlayout.setStretchFactor(self.pathwaylist, 9)
        self.pathwaycount = QLabel('Count: ')
        self.pathwaylistlayout.addWidget(self.pathwaycount)
        self.pathwaylistlayout.setStretchFactor(self.pathwaycount, 1)
        self.main_layout.addLayout(self.pathwaylistlayout)
        self.main_layout.setStretchFactor(self.pathwaylistlayout, 4)

        self.processlist = QListWidget()
        #self.main_layout.addWidget(self.processlist)
        self.processcount = QLabel('Count: ')
        #self.main_layout.addWidget(self.processcount)

    def dropdownmenu(self):
        self.spinsystem_menu = QComboBox()
        self.spinsystem_menu.addItems(list(self.spinsystems_widget.spinsystems.keys()))
        self.left_layout.addWidget(self.spinsystem_menu)
        def on_currentIndexChanged():
            self.spinsystems_widget.built_motifs = {}
            self.builtmotiflist.clear()
            self.builtmotifcount.setText('Count: ')
            self.pathwaylist.clear()
            self.pathwaycount.setText('Count: ')
            self.spinsystems_widget.displayspinsystem(self.spinsystem_menu.currentText())
        self.spinsystem_menu.currentIndexChanged.connect(on_currentIndexChanged)

    def structureoptions(self):
        self.structure_menu = {'layout': QVBoxLayout()}
        text = {'metabolites': 'Metabolites', 'motif_0': 'Motif_0', 'motif_1': 'Motif_1', 'motif_2': 'Motif_2', 'submotif_1': 'Sub-Motif_1', 'submotif_2': 'Sub-Motif_2', 'submotif_3': 'Sub-Motif_3', 'submotif_4': 'Sub-Motif_4'}
        for opt in ['metabolites', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4']:
            self.structure_menu[opt] = QCheckBox(text[opt])
            self.structure_menu[opt].setChecked(True)
            self.structure_menu['layout'].addWidget(self.structure_menu[opt])
        self.left_layout.addLayout(self.structure_menu['layout'])

    def speciesmenu(self):
        self.species_menu = QComboBox()
        specieslist = self.spinsystems_widget.pathwaymapper.Reactome.specieslist
        specieslist.remove('Homo sapiens')
        specieslist = ['Homo sapiens']+sorted(specieslist, key=lambda x: x)
        self.species_menu.addItems(specieslist)
        self.left_layout.addWidget(self.species_menu)
    
    # forward our widget events
    #def keyPressEvent(self, event):         
        #self.main_tabs.currentWidget().keyPressEvent(event)

    #def mousePressEvent(self, event):
        #self.main_tabs.currentWidget().mousePressEvent(event)

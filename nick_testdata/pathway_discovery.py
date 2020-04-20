import pandas, numpy

with open('nick_testdata.xlsx', 'rb') as xlsx:
    data = pandas.read_excel(xlsx)

data_set = numpy.array([data['NT_CFvsWT'], data['Unnamed: 1'], data['Unnamed: 2']])

import csv

names_to_colmarid = {}
with open('../updating_Reactome_molecules/name_to_colmarid.csv') as csvin:
    csvin.readline()
    for row in csv.reader(csvin):
        names_to_colmarid[row[0]] = row[1]

metabolites = {}
for i in range(len(data_set[0])):
    metabolite, fold_change = data_set[0][i], data_set[2][i]
    if metabolite in names_to_colmarid:
        colmar = names_to_colmarid[metabolite]
        metabolites[metabolite] = {'COLMAR': colmar, 'Fold Change': fold_change}

import sys

sys.path.append('../scripts')
from neo4j_connect import ReactomeServer
Reactome = ReactomeServer()

with Reactome._driver.session() as db:
    names = db.run('MATCH (m:ReferenceMolecule) RETURN split(m.displayName, " ")').value()
    for n in names:
        if 'L-lysinium(1+)' in n[0]:
            print(n)

    for metabolite in metabolites:
        colmar = metabolites[metabolite]['COLMAR']
        pathways = db.run('MATCH (m:ReferenceMolecule) WHERE "'+colmar+'" IN m.COLMAR RETURN m.displayName').value()
        print(metabolite, colmar, pathways)
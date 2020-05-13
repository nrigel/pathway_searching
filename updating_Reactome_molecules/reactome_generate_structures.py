# parts 1 and 3 require py 3.6, part 2 requires 2.7

def part1():

    import sys, openbabel, pybel, pickle, csv

    sys.path.append('../scripts')
    from neo4j_connect import ReactomeServer
    Reactome = ReactomeServer()

    reactome_to_chebi = {}
    with Reactome._driver.session() as db:
        reactome_compounds = db.run('MATCH (m:ReferenceMolecule) RETURN m.displayName').value()

        for compound in reactome_compounds:
            chebi = ''
            for i in range(compound.find('ChEBI'), len(compound)):
                if compound[i] == ']' or compound[i] == ' ':
                    break
                chebi = chebi+compound[i]
            chebi = 'CHEBI:'+chebi.split(':')[1]

            reactome_to_chebi[compound] = chebi

    with open('reactome_to_chebi.csv', 'w') as pkl:
        writer = csv.writer(pkl)
        for k in reactome_to_chebi:
            writer.writerow([k, reactome_to_chebi[k]])
        

def part2():
    import sys, pickle, openbabel, pybel, csv
    sys.path.append('../scripts')
    from svg_drawer import MOLDrawer

    D = MOLDrawer()
    reactome_to_chebi = {}
    with open('reactome_to_chebi.csv', 'r') as pkl:
        for row in csv.reader(pkl):
            reactome_to_chebi[row[0]] = row[1]

    chebi_to_reactome = {reactome_to_chebi[k]: k for k in reactome_to_chebi}

    reactome_to_svg = {}
    for chebi in pybel.readfile("sdf", "chebi_hmdb_conversion/ChEBI_complete.sdf"):
        if chebi.data['ChEBI ID'] in chebi_to_reactome:
            print(chebi.data['ChEBI ID'])
            mol = chebi.write('mol')
            try:
                reactome_to_svg[chebi_to_reactome[chebi.data['ChEBI ID']]] = D.Draw(mol)
            except Exception as e:
                print(e)

    with open('reactome_SVGs.pkl', 'wb') as pkl:
        pickle.dump(reactome_to_svg, pkl, pickle.HIGHEST_PROTOCOL)

def part3():
    import sys, pickle
    with open('reactome_SVGs.pkl', 'rb') as pkl: 
        reactome_to_svg = pickle.load(pkl)

    sys.path.append('../scripts')
    from neo4j_connect import ReactomeServer
    Reactome = ReactomeServer()

    with Reactome._driver.session() as db:
        for m in reactome_to_svg:
            db.run('MATCH (m:ReferenceMolecule) WHERE m.displayName = "'+m+'" SET m.SVG = "'+reactome_to_svg[m]+'"')

part3()
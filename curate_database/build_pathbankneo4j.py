
# node csvs: 'dbId:ID',descriptors,':LABEL'
# edge csvs: ':START_ID',descriptors,':END_ID',':TYPE'
neo4j_admin_path = '/Applications/neo4j-community-3.5.14_pathbank/bin/'
path = 'pathbank_data/'
pathbank_nodes = {'Pathway': {'CSV': 'pathbank_pathways.csv', 'id_index': 0, 'species_index': 4, 'header': ['SMPDB_ID:ID', 'PW_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Description', ':LABEL']},
                'Metabolite': {'CSV': 'pathbank_all_metabolites.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'PW_ID:ID', 'Metabolite_Name', 'HMDB_ID', 'KEGG_ID', 'ChEBI_ID', 'DrugBank_ID', 'CAS', 'Formula', 'IUPAC', 'SMILES', 'InChI', 'InChIKey', ':LABEL']},
                'Protein': {'CSV': 'pathbank_all_proteins.csv', 'id_index': 4, 'species_index': 3, 'header': ['SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Uniprot_ID:ID', 'Protein_Name', 'HMDBP_ID', 'DrugBank_ID', 'GenBank_ID', 'Gene_Name', 'Locus', ':LABEL']}}

node_data, node_csvs, edge_csvs = {}, [], []
pathway_species = {}
import csv, json

for node in ['Metabolite', 'Protein', 'Pathway']:
    print('Modifying '+node+' nodes...')
    node_data[node] = {}
    node_csvs.extend([path+'neo4jnodes_'+node+'.csv'])
    desc = {}

    with open(path+pathbank_nodes[node]['CSV']) as csvin, open(node_csvs[-1], 'w') as csvout:
        writer = csv.writer(csvout)
        writer.writerow(pathbank_nodes[node]['header'])
        csvin.readline()
        for row in csv.reader(csvin):
            ID = row[pathbank_nodes[node]['id_index']]
            write_row = row+[node]

            if node != 'Pathway':
                pathway_species[row[0]] = row[pathbank_nodes[node]['species_index']]
            else:
                if ID not in pathway_species:
                    species = 'Homo sapiens'
                else:
                    write_row.insert(pathbank_nodes[node]['species_index'], pathway_species[ID])

            write_row = [json.dumps(r.replace('"', '')) for r in write_row]
            write_row = ','.join(write_row)+'\n'
            if node == 'Pathway':
                csvout.write(write_row)
            else:
                if node == 'Protein':
                    d = tuple(row[4:11])
                    if ID in node_data[node]:
                        m = len([i for i in d if i == ''])
                        M = len([i for i in desc[ID] if i == ''])
                        if m < M:
                            desc[ID] = d
                            node_data[node][ID] = write_row
                    else:
                        desc[ID] = d
                        node_data[node][ID] = write_row
                else:
                    node_data[node][ID] = write_row
        if node != 'Pathway':
            for ID in node_data[node]:
                csvout.write(node_data[node][ID])


# add in details like Reactome...
    # already have HMDB, CHEBI IDs
    # link SVG structures, COLMAR IDs, anything else we want, do that here
# What we need for a database...
    # Unified nodes: Pathways, Metabolites, Proteins
    # Relationships: Pathways -> stepConversion -> BioChemicalReaction -> participantStoichiometry -> Stoichiometry -> physicalEntity (protein or metabolite ref), stoichiometricCoefficient
    #                                                                  -> displayName, eft, right, name, spontaneous?, stoichiometry, conversionDirection

neo4j_admin_command = [neo4j_admin_path+'neo4j-admin import ']
for file in node_csvs:
    neo4j_admin_command.extend(['--nodes='+file])
for file in edge_csvs:
    neo4j_admin_command.extend(['--relationships='+file])
print(' '.join(neo4j_admin_command))


# /Applications/neo4j-community-3.5.14_pathbank/bin/neo4j-admin import  --nodes=pathbank_data/neo4jnodes_Pathway.csv --nodes=pathbank_data/neo4jnodes_Metabolite.csv --nodes=pathbank_data/neo4jnodes_Metabolite_aliases.csv --nodes=pathbank_data/neo4jnodes_Protein.csv --nodes=pathbank_data/neo4jnodes_Protein_aliases.csv --relationships=pathbank_data/neo4jedges_Metabolite_aliases.csv --relationships=pathbank_data/neo4jedges_Protein_aliases.csv

# node csvs: 'dbId:ID',descriptors,':LABEL'
# edge csvs: ':START_ID',descriptors,':END_ID',':TYPE'
neo4j_admin_path = '/Applications/neo4j-community-3.5.14_pathbank/bin/'
path = 'pathbank_data/'
pathbank_nodes = {'Pathway': {'CSV': 'pathbank_pathways.csv', 'id_index': 0, 'owl_index': 0, 'header': ['dbId:ID', 'SMPDB_ID', 'PW_ID', 'Pathway_Name', 'Pathway_Subject', 'Description', ':LABEL']},
                'Metabolite': {'CSV': 'pathbank_all_metabolites.csv', 'id_index': 4, 'owl_index': 0, 'header': ['dbId:ID', 'SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'PW_ID', 'Metabolite_Name', 'HMDB_ID', 'KEGG_ID', 'ChEBI_ID', 'DrugBank_ID', 'CAS', 'Formula', 'IUPAC', 'SMILES', 'InChI', 'InChIKey', ':LABEL']},
                'Protein': {'CSV': 'pathbank_all_proteins.csv', 'id_index': 4, 'owl_index': 0, 'header': ['dbId:ID', 'SMPDB_ID', 'Pathway_Name', 'Pathway_Subject', 'Species', 'Uniprot_ID', 'Protein_Name', 'HMDBP_ID', 'DrugBank_ID', 'GenBank_ID', 'Gene_Name', 'Locus', ':LABEL']}}

IDs, OWLs, node_csvs, edge_csvs = {}, {}, [], []

import csv, pickle, json

for node in pathbank_nodes:
    print('Modifying '+node+' nodes...')
    IDs[node], OWLs[node] = {}, {}
    node_csvs.extend([path+'neo4jnodes_'+node+'.csv'])
    with open(path+pathbank_nodes[node]['CSV']) as csvin, open(node_csvs[-1], 'w') as csvout:
        writer = csv.writer(csvout)
        writer.writerow(pathbank_nodes[node]['header'])
        line_counter = -1
        for line in csvin:
            line_counter += 1
        csvin.seek(0)
        csvin.readline()
        
        row_counter = 0
        for row in csv.reader(csvin):
            row_counter += 1
            dbId = node+'_'+str(row_counter).zfill(len(str(line_counter))+1)
            ID = row[pathbank_nodes[node]['id_index']]
            OWL_ID = row[pathbank_nodes[node]['owl_index']]+'-'+ID
            if ID not in IDs[node]:
                IDs[node][ID] = []
            if OWL_ID not in OWLs[node]:
                OWLs[node][OWL_ID] = []
            IDs[node][ID].extend([dbId])
            OWLs[node][OWL_ID].extend([dbId])

            write_row = [dbId]+row+[node]
            write_row = [json.dumps(r.replace('"', '')) for r in write_row]
            write_row = ','.join(write_row)+'\n'
            csvout.write(write_row)

    IDs[node] = {i: IDs[node][i] for i in IDs[node] if len(IDs[node][i]) > 1}
    if len(IDs[node]) == 0:
        continue
    
    print('\tRelating replicate nodes ('+str(len(IDs[node]))+')...')

    node_csvs.extend([path+'neo4jnodes_'+node+'_aliases.csv'])
    edge_csvs.extend([path+'neo4jedges_'+node+'_aliases.csv'])
    with open(node_csvs[-1], 'w') as nodeout, open(edge_csvs[-1], 'w') as edgeout:
        nwriter, ewriter = csv.writer(nodeout), csv.writer(edgeout)
        nwriter.writerow(['dbId:ID', pathbank_nodes[node]['header'][pathbank_nodes[node]['id_index']+1], ':LABEL'])
        ewriter.writerow([':START_ID', ':END_ID', ':TYPE'])

        for ID in IDs[node]:
            row_counter += 1
            dbId = node+'_'+str(row_counter).zfill(len(str(line_counter))+1)
            nwriter.writerow([dbId, ID, node+'Alias'])
            for child in IDs[node][ID]:
                ewriter.writerow([child, dbId, 'nodeAlias'])

with open('OWL_IDs.pkl', 'wb') as pkl:
    pickle.dump(OWLs, pkl, pickle.HIGHEST_PROTOCOL)


# OWL_IDs.pkl: {Node: {PathwayID-NodeID: [dbIds]}}
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
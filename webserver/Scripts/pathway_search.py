#!/Users/nick/Documents/GitHub/motif_builder/py35env/bin/python
#!/opt/anaconda2/envs/P3/bin/python
#!/usr/bin/env python3
# Python 3.6.10 :: Anaconda, Inc.#!/opt/anaconda2/envs/P3/bin/python

print("Content-type:text/html\r\n\r\n") # print HTML header

# Set up Neo4j connections
# set parameters for Bolt port, user, and password
node_format = [{'size': 30, 'color': ['#E50002', '#E02C00', '#DC5900', '#D88500', '#D4AF00', '#C8CF00', '#99CB00', '#6CC700', '#41C300', '#18BF00'], 'shape': 'ellipse'}, 
                {'size': 30, 'color': '#0074BF', 'shape': 'square'}]


from neo4j_connect import Neo4jServer

NMRT, PathBank = Neo4jServer('NMRT'), Neo4jServer('PathBank')

from session_save import SaveSession
import cgi, cgitb, math, json, csv, traceback, time, io
start = time.time()
# Step 1: Organize data from HTML form
form = cgi.FieldStorage() # Create instance of FieldStorage 

# Get data from fields
cliques = []
for row in csv.reader([form.getvalue(c) for c in form.keys() if c[:4] == 'cbox']):
    ROW = [row[0]]+[r[1:] for r in row[1:]]
    cliques.extend([tuple(ROW)])

prelim = False
if form.getvalue("search_opt") == 'true':
    prelim = True

saveinfo = json.loads(form.getvalue('saveinfo'))
cached_results = json.loads(form.getvalue('cached_results'))
inputOptions = cached_results['inputOptions']

structure_opts = [o for o in ['metabolites', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4'] if form.getvalue(o)]
species = form.getvalue('specieslist')
count_cutoff = form.getvalue('cutoff_count') # change to dict of cutoffs for all structure opts
clique_cutoff = form.getvalue('clique_cutoff')
# Read in data file
if not cached_results['filedata']:
    datafile = form['userfile'].file
    data_set = [line.rstrip().decode('utf-8') for line in datafile.readlines()]
else:
    data_set = cached_results['filedata'].split('\n')

inputted_file = {'data_set': data_set, 'filename': form['userfile'].filename}
data_set = [row for row in csv.reader(data_set)]

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

try:
    # Step 2: Match inputted metabolites to NMRT DB to get (sub-)motif info and COLMAR IDs (if needed)
    metabolites = {}
    fold_changes, p_values = [], []
    #start = time.time()
    with NMRT._driver.session() as db:
        for i in range(1, len(data_set)):
            metabolite = data_set[i][0]
            
            if isfloat(data_set[i][1]):
                p_value = float(data_set[i][1])
            else:
                p_value = 1
            if isfloat(data_set[i][2]):
                fold_change = float(data_set[i][2])
            else:
                fold_change = 1

            fold_changes.extend([fold_change])
            p_values.extend([p_value])

            colmar_matches = db.run('MATCH (m:metabolites) WHERE m.COLMARm = "'+metabolite+'" OR m.COLMARm = "'+metabolite+'_1" OR m.COLMARm = "'+metabolite+'_2" OR m.Name = "'+metabolite+'" OR m.isoName = "'+metabolite+'" RETURN m.COLMARm').value()
            metabolites[metabolite] = {'metabolites': colmar_matches, 'Fold Change': fold_change, 'P-Value': p_value, 'Metabolite': metabolite}
            
            for opt in [o for o in structure_opts if o != 'metabolites']:
                metabolites[metabolite][opt] = set()
                for match in colmar_matches:
                    for s in list(set(db.run('MATCH (m:metabolites)-[*..2]-(M:'+opt+') WHERE m.COLMARm = "'+match+'" RETURN M.SMILES').value())):
                        metabolites[metabolite][opt].add(s)

    # Step 3: Search for pathways that contain matching metabolites, motifs, and/or sub-motifs
    pathway_search = {'metabolites': {}, 'pathways': {}, 'pathways_to_return': {}, 'pathway_graphs': {}}
    pathway_results = {'pathway_data': {}, 'pathway_list': {}, 'pathway_cliques': {}, 'structure_options': structure_opts}
    
    with PathBank._driver.session() as db:
        pathway_species = {p[0]: p[1] for p in db.run('MATCH (p:Pathway) RETURN ([p.SMPDB_ID, p.Species])').value()}
        # iterate thru all options given in the form
        for opt in structure_opts:
            pathwaylist = [] # save pathway names in order of most to least matching metabolites

            pathway_search['metabolites'][opt] = {} # reactome metabolite to colmar metabolite
            pathway_search['pathways'][opt] = {} # pathways to set of metabolites in them
            pathway_search['pathways_to_return'][opt] = [] # pathways that are kept based on the count cutoff
            pathway_search['pathway_graphs'][opt] = {} # pathway to nx graphs that will contain the metabolites and reaction connectivity
            pathway_results['pathway_data'][opt] = {}
            pathway_results['pathway_cliques'][opt] = {'keys': [], 'values': [], 'clique_pathways': {}}

            for metabolite in metabolites: # iterate thru user inputted metabolites
                for s in metabolites[metabolite][opt]: # iterate thru structure options
                    if opt == 'metabolites':
                        key = 'COLMARm'
                    else:
                        key = opt
                    for (m, plist) in db.run('MATCH (m:COLMARm) WHERE "'+s+'" IN m.'+key+' RETURN ([m.PW_ID, m.SMPDB_ID])').value():
                        for p in plist:
                            if p not in pathway_species or pathway_species[p] != species:
                                continue
                            # m = reactome metabolite, r = reaction, p = pathway
                            # keep m to metabolite alias
                            if m not in pathway_search['metabolites'][opt]:
                                pathway_search['metabolites'][opt][m] = set()
                            pathway_search['metabolites'][opt][m].add(metabolite) # metabolite matches are returned later
                            # add metabolite to pathway dict set to keep track of occurences
                            if p not in pathway_search['pathways'][opt]:
                                pathway_search['pathways'][opt][p] = set() # metabolite list
                            #pathway_search['pathways'][opt][p].add(metabolite)
                            pathway_search['pathways'][opt][p].add(s)

            # get pathways to return as result
            for pathway in pathway_search['pathways'][opt]:
                metcount = len(pathway_search['pathways'][opt][pathway])
                if metcount >= int(count_cutoff):
                    mets = tuple([opt]+sorted([m for m in pathway_search['pathways'][opt][pathway]], key=lambda x: x))
                    if mets not in pathway_results['pathway_cliques'][opt]['clique_pathways']:
                        pathway_results['pathway_cliques'][opt]['clique_pathways'][mets] = set()
                    pathway_results['pathway_cliques'][opt]['clique_pathways'][mets].add(pathway)

            for pathway in pathway_search['pathways'][opt]:    
                metcount = len(pathway_search['pathways'][opt][pathway])
                if metcount >= int(count_cutoff):
                    mets = tuple([opt]+sorted([m for m in pathway_search['pathways'][opt][pathway]], key=lambda x: x))
                    pathcount = len(pathway_results['pathway_cliques'][opt]['clique_pathways'][mets])
                    if mets not in pathway_results['pathway_cliques'][opt]['keys']:
                        pathway_results['pathway_cliques'][opt]['keys'].extend([mets])
                        pathway_results['pathway_cliques'][opt]['values'].extend([pathcount])
                    if len(cliques) == 0 and pathcount > int(clique_cutoff): # cutoff given in form (need to change to dict instead)
                        continue
                    if len(cliques) > 0 and mets not in cliques:
                        continue
                    pathway_search['pathways_to_return'][opt].extend([pathway])
            del pathway_results['pathway_cliques'][opt]['clique_pathways'] # have to do this because you cannot JSON dump a dict with tuple keys
            
            #print('metabolite matching finished', time.time()-start)
            
            if prelim:
                continue
            
            reaction_labels = ['BiochemicalReaction', 'Interaction', 'MolecularInteraction', 'Transport', 'TransportWithBiochemicalReaction']

            # iterate thru pathways that meet the cutoffs
            for pathway in pathway_search['pathways_to_return'][opt]:
                pathway_search['pathway_graphs'][opt][pathway] = [{}, [], {}, []] # nodedict, nodes, edges
                G = pathway_search['pathway_graphs'][opt][pathway]

                result = db.run('MATCH (m:Metabolite)-[er]-(r:Reaction)-[he:hasEvent]-(p:Pathway {SMPDB_ID: "'+pathway+'"}) RETURN ([[x in keys(m) WHERE not x in ["SMPDB_ID"]| [x, m[x] ] ], type(er), r, he, p])').value()
                reactions = {} # reaction to set of metabolite, input/output pairs
                pathway_name = pathway
                for (M, er, r, he, p) in result: # iterate thru results
                    m = {i[0]: i[1] for i in M}
                    if p['Pathway_Name']:
                        pathway_name = p['Pathway_Name'] # save pathway stuff as the name instead of the id
                    if 'PW_ID' in m:
                        m_ID = m['PW_ID']
                    else:
                        m_ID = m['Biopax_ID']
                    if m_ID not in G[0]: # add node to G with all of the Reactome info
                        node_dict = dict(m) # convert reactome node object to dict
                        node_dict['id'] = m_ID # need this for cy plot or it will make up its own
                        node_dict['type'] = 'Metabolite' # need this for cy options later in the js
                        if 'Metabolite_Name' in m:
                            node_dict['displayName'] = m['Metabolite_Name']

                        if 'COLMARm' not in node_dict:
                            node_dict['COLMARm'] = []
                        
                        if m_ID in pathway_search['metabolites'][opt]: # means that we can add in user submitted data, such as fold change, p-value 
                            # figure this out later
                            colmar = list(pathway_search['metabolites'][opt][m_ID])[0]
                            node_dict['Metabolite Matches'] = list(pathway_search['metabolites'][opt][m_ID]) # save this to show later
                            node_dict['Fold Change'], node_dict['P-Value'] = metabolites[colmar]['Fold Change'], metabolites[colmar]['P-Value']
                            # fold change determines color selected from gradient
                            if max(fold_changes) != min(fold_changes):
                                fcfactor = (node_dict['Fold Change']-min(fold_changes))/(max(fold_changes)-min(fold_changes))*10
                                color = 0
                                for i in range(0, 10):
                                    if abs(fcfactor-i) < abs(fcfactor-color):
                                        color = i
                                node_dict['color'] = node_format[0]['color'][color]
                            else:
                                node_dict['color'] = node_format[1]['color']
                            # p-value determines size scaling
                            if node_dict['P-Value'] > 0:
                                node_dict['size'] = abs(math.log(node_dict['P-Value']))*node_format[0]['size']
                            else:
                                node_dict['size'] = node_format[0]['size']
                            
                        else:
                            node_dict['size'] = node_format[1]['size']
                            node_dict['color'] = node_format[1]['color']

                        if 'COLMARm' not in node_dict or len(node_dict['COLMARm']) == 0:
                            node_dict['shape'] = node_format[1]['shape']
                        else:
                            node_dict['shape'] = node_format[0]['shape']
                        
                        G[0][m_ID] = node_dict
                        G[1].extend([{'data': node_dict, 'group': 'nodes'}])

                    if r['Biopax_ID'] not in reactions: # keep reaction info
                        rdict = dict(r)
                        rdict['displayName'] = r['displayName']
                        rdict['stoichiometry'], rdict['order'] = he['stoichiometry'], he['order']
                        reactions[r['Biopax_ID']] = {'metabolites': set(), 'node': rdict} # save node object to return to js
                    reactions[r['Biopax_ID']]['metabolites'].add((m_ID, er))

                for r in reactions: # iterate thru reactions
                    for (m, er) in reactions[r]['metabolites']:
                        for (M, ER) in reactions[r]['metabolites']:
                            if m == M or er == ER: # don't want edges for two inputs or outputs, or for a node to itself
                                continue
                            edge_dict = {'data': {'reaction': reactions[r]['node']}, 'group': 'edges'}
                            if er == 'input': # directed edges
                                edge_dict['data']['source'] = m
                                edge_dict['data']['target'] = M
                            if er == 'output':
                                edge_dict['data']['source'] = M
                                edge_dict['data']['target'] = m
                            if (edge_dict['data']['source'], edge_dict['data']['target'], edge_dict['data']['reaction']['Biopax_ID']) not in G[2]:
                                G[2][(edge_dict['data']['source'], edge_dict['data']['target'], edge_dict['data']['reaction']['Biopax_ID'])] = edge_dict
                                G[3].extend([edge_dict])

                # include metabolite count in pathway name for our pathway list in the html
                pathway_name = pathway_name+' ('+str(len(pathway_search['pathways'][opt][pathway]))+')'
                # add to this list so that we can order it
                pathwaylist.extend([(pathway_name, len(pathway_search['pathways'][opt][pathway]))])
                # save to dict for js to read in
                pathway_results['pathway_data'][opt][pathway_name] = {'nodes': G[1], 'edges': G[3]}

            # save pathway list to results dict
            pathway_results['pathway_list'][opt] = [p[0] for p in sorted(pathwaylist, key=lambda x: (-x[1], x[0]))]

        for opt in pathway_results['structure_options']:
            indices = []
            for clique in sorted(pathway_results['pathway_cliques'][opt]['keys'], key=lambda x: -pathway_results['pathway_cliques'][opt]['values'][pathway_results['pathway_cliques'][opt]['keys'].index(x)]):
                n = pathway_results['pathway_cliques'][opt]['keys'].index(clique)
                indices.extend([n])
            pathway_results['pathway_cliques'][opt]['keys'] = [pathway_results['pathway_cliques'][opt]['keys'][i] for i in indices]
            pathway_results['pathway_cliques'][opt]['values'] = [pathway_results['pathway_cliques'][opt]['values'][i] for i in indices]
        
        pathway_results['inputted_file'] = inputted_file # save user file for session loading later
        pathway_results['inputOptions'] = cached_results['inputOptions'] # include our user input options in our save
        pathway_results['saveinfo'] = SaveSession(db, pathway_results, saveinfo) # save here instead of making new request
        


    print(json.dumps(pathway_results)) # convert to JSON and print results back to js
    
    PathBank.close()
    NMRT.close()
except Exception as e:
    print(traceback.format_exc())
    exit()

### To Do List ###
# try collapsible lists again - this time save list contents to var, then clear list when button is pressed
# write function to remove white space from svg
# figure out why isolated nodes in a pathway are thrown far away and correct that
# run examples from literature
# link unique pathways to infection type (biofilm vs planktonic)


# Can we scrape PathBank data and curate our own database???

# Work on saving session, results


### Rafael suggestions that are not yet possible to implement ###
# proteomics, RNAseq -> ReferenceRNASequence
# missing metabolites in sample cohorts (COLMARq; p-values)

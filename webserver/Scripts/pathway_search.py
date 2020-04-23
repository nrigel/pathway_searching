#!/opt/anaconda2/envs/P3/bin/python

import cgi, cgitb 

# Create instance of FieldStorage 
form = cgi.FieldStorage() 

print("Content-type:text/html\r\n\r\n")
print('<html>')
print('<head>')
print('</head>')
print('<body>')

# Get data from fields
#print(form.keys())
datafile = form.getvalue('userfile')
structure_opts = {o: form.getvalue(o) for o in ['metabolites', 'motif_0', 'motif_1', 'motif_2', 'submotif_1', 'submotif_2', 'submotif_3', 'submotif_4']}
species = form.getvalue('specieslist')
count_cutoff = form.getvalue('cutoff_count')
#print(datafile)
print(structure_opts)
print(species)
print(count_cutoff)


# Need to figure out how to keep selected file on page reload
# Need to call search script from here with params
# Need to return search results
# Need to figure out how I want to display results


print('</body>')
print('</html>')
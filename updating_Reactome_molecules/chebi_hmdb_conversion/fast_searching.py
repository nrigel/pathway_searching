import openbabel, pybel, pickle

#obConversion.SetInAndOutFormats('fs', 'sdf')
#index = openbabel.OBMol()
#obConversion.ReadString(ob_mol, "C=CC#N")

conv = openbabel.OBConversion()

#conv.OpenInAndOutFiles("HMDB_structures.sdf","HMDB_index.fs")
#conv.SetInAndOutFormats("sdf","fs")
#conv.Convert()

chebi_to_hmdb = {}
for chebi in pybel.readfile("sdf", "ChEBI_complete.sdf"):
    chebi_to_hmdb[chebi.data['ChEBI ID']] = []
    
    conv.OpenInAndOutFiles("HMDB_index.fs", "results.sdf")
    conv.SetInAndOutFormats("fs", "sdf")
    conv.AddOption("s", conv.GENOPTIONS, chebi.write('can'))
    conv.AddOption("e", conv.INOPTIONS)
    conv.Convert()
    for hmdb in pybel.readfile("sdf", "results.sdf"):
        chebi_to_hmdb[chebi.data['ChEBI ID']].extend([hmdb.data['HMDB_ID']])

with open('CHEBI_TO_HMDB.pkl', 'wb') as pkl:
    pickle.dump(chebi_to_hmdb, pkl, pickle.HIGHEST_PROTOCOL)

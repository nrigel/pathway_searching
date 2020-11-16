
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

def Draw(SMILES, annotateAtoms={}, molSize=(450,150)):
    
    if annotateAtoms:
        for key in list(annotateAtoms.keys()):
            if type(key) != int:
                annotateAtoms[int(key)] = annotateAtoms[key]
                del annotateAtoms[key]

    mol = Chem.MolFromSmiles(SMILES)
    if not mol:
        return
    mc = Chem.Mol(mol.ToBinary())
    Chem.rdDepictor.Compute2DCoords(mc)
    
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    opts = drawer.drawOptions()
    for atom in annotateAtoms:
        opts.atomLabels[atom-1] = annotateAtoms[atom]
    
    bond_colors = {bond.GetIdx(): (0.91,0.4,0.4) for bond in mc.GetBonds() if bond.GetBeginAtomIdx()+1 in annotateAtoms and bond.GetEndAtomIdx()+1 in annotateAtoms}

    drawer.DrawMolecule(mc,highlightAtoms=[atom-1 for atom in annotateAtoms.keys()],highlightAtomColors={atom-1: (0.91,0.4,0.4) for atom in annotateAtoms}, highlightBonds=bond_colors)
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

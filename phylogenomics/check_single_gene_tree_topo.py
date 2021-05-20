import glob
import ete3
import os 
    
def get_mono_clades(tree, taxon):
    seeds = set([l for l in tree.get_leaves() if l.clade == taxon])
    nodes = set()
    for s in seeds:
        n = s
        while all([l.clade == taxon for l in n.up.get_leaves()]):
            n = n.up
        nodes.add(n)
    return nodes
    
def get_sister_id(node):
    if node.up.get_children()[0] == node:
        return set([l.clade for l in node.up.get_children()[1]])
    elif node.up.get_children()[1] == node:
        return set([l.clade for l in node.up.get_children()[0]])


with open('sisters_picozoa_sgt.csv', 'w') as out:
    for f in glob.glob("trees_for_fabien_renamed/*.new"):
        clst = os.path.basename(f).replace('.new', '')
        tree = ete3.PhyloTree(f, format=2)
        for l in tree.iter_leaves():
            l.add_feature(pr_name='clade', pr_value=l.name.replace("'", "").split('_')[0])
        clades = get_mono_clades(tree, 'Picozoa')
        for i, n in enumerate(clades, 1):
            sister = get_sister_id(n)
            if sister:
                print("{}\tclade {}\t{}".format(clst, i, ";".join(list(sister))), file=out)

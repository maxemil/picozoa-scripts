import ete3
import glob


def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor

def get_clade(tree, clades):
    prr = [l for l in tree.get_leaves() if any([c in l.name for c in clades])]
    if prr:
        return get_ancestor(prr)


for f in glob.glob('sgt_bootstraps/*.treefile'):
    tree = ete3.PhyloTree(f, format=0)
    pico = get_clade(tree, ['Picozoa'])
    prr = ['Picozoa', 'Rhodelphis', 'Rhodophyta']
    if pico and not pico.is_root():
        anc = pico.up
        leaf_names = [l.name for l in anc.get_leaves()]
        while all([any([c in l for c in prr]) for l in leaf_names]):
            print(f, set([l.split('_')[0] for l in leaf_names]))
            anc = anc.up
            leaf_names = [l.name for l in anc.get_leaves()]

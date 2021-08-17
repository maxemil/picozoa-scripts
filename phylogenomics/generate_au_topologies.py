import ete3

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor

def get_clade(tree, clades):
    prr = [l for l in tree.get_leaves() if any([c in l.name for c in clades])]
    return get_ancestor(prr)

def make_sisters(tree, clade1, clade2):
    anc = get_clade(tree, clade1)
    anc_up = anc.up
    anc.detach()
    anc_up.delete()

    sis = get_clade(tree, clade2)
    sis_up = sis.up
    sis.detach()
    
    new_node = ete3.PhyloNode()
    sis_up.add_child(new_node)
    new_node.add_child(sis)
    new_node.add_child(anc)


topologies = [(['Picozoa'], ['Rhodelphis', 'Rhodophyta']),
(['Picozoa'], ['Rhodophyta']),
(['Picozoa'], ['Rhodelphis']),
(['Picozoa'], ['Viridiplantae', 'Glaucophyta']),
(['Picozoa'], ['Glaucophyta']),
(['Picozoa'], ['Viridiplantae']),
(['Picozoa'], ['Archaeplastida']),
(['Picozoa'], ['Telonemia']),
(['Picozoa'], ['Telonemia', 'Rhizaria', 'Stramenopila', 'Alveolata']),
(['Picozoa'], ['Cryptista']),
(['Picozoa', 'Cryptista'], ['Rhodophyta', 'Rhodelphis']),
(['Picozoa', 'Cryptista'], ['Rhodophyta']),
(['Picozoa', 'Cryptista'], ['Viridiplantae', 'Glaucophyta']),
(['Picozoa', 'Cryptista'], ['Glaucophyta']),
(['Picozoa', 'Cryptista'], ['Viridiplantae'])]

tree = ete3.PhyloTree("orig_topology.new", format=0)

with open('all_topologies2test.new', 'w') as out:
    for c1, c2 in topologies:
        # print("({}),({})".format(','.join([c for c in c1]), ','.join([c for c in c2])))
        make_sisters(tree, c1, c2)
        print(tree.write(format=9), file=out)

import ete3

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor


tree = ete3.PhyloTree('astral/astral_tree.main_tree', quoted_node_names=False)
treepp = ete3.PhyloTree('astral/astral_tree_pp.new', quoted_node_names=False)


counter = 0
c2supp = {}
for n in tree.traverse():
    if not n.is_leaf():
        kids = [c.name for c in n.get_leaves()]
        pp_n = get_ancestor([treepp.get_leaves_by_name(k)[0] for k in kids])
        c2supp[counter] = "{}/{}".format(round(n.support, 1), pp_n.support)
        n.support = counter
        counter += 1


tree_str = tree.write()
for c,s in c2supp.items():
    tree_str = tree_str.replace(")" + str(c) + ":", ")'" + s + "':")


with open('astral/astral_tree_combined.new', 'w') as out:
    print(tree_str, file=out)

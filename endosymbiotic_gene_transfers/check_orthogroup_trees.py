import glob
import ete3
from collections import defaultdict
import shutil
import os

plastid_groups = ['Apicomplexa',
                 'Chlorarachniophyceae',
                 'Chloroplastida',
                 'Colpodellida',
                 'Cryptophyceae',
                 'Cyanobacteria',
                 'Dinoflagellata',
                 'Euglenida',
                 'Glaucophyta',
                 'Haptophyta',
                 'Ochrophyta',
                 'Paulinella',
                 'Picozoa',
                 'Rhodelphis',
                 'Rhodophyta']


def find_sisters(tree, mincount, groupOI):
    for n in tree.traverse():
        if not n.is_leaf():
            taxa_count_c1 = count_taxa(n.get_children()[0])
            plastid_count_c1 = count_plastid(taxa_count_c1)
            taxa_count_c2 = count_taxa(n.get_children()[1])
            plastid_count_c2 = count_plastid(taxa_count_c2)
            if is_clade_monophyletic(plastid_count_c1, 'plastid') and taxa_count_c1[groupOI] >= mincount and is_clade_monophyletic(taxa_count_c2, 'Cyanobacteria') or \
                 is_clade_monophyletic(plastid_count_c2, 'plastid') and taxa_count_c2[groupOI] >= mincount and is_clade_monophyletic(taxa_count_c1, 'Cyanobacteria'):
                 return True
    return False

def count_plastid(taxa_count):
    plastid_count = {'plastid':0, 'other':0}
    for k,v in taxa_count.items():
        if k in plastid_groups:
            plastid_count['plastid'] += v
        else:
            plastid_count['other'] += v
    return plastid_count

def is_clade_monophyletic(taxa_count, group, impurity=0.1):
    if taxa_count[group] >= (1-impurity) * sum(taxa_count.values()):
        return True
    else:
        return False

def count_taxa(node):
    taxa_count = defaultdict(int)
    for l in node.iter_leaves():
        taxa_count[l.clade] += 1
    return taxa_count
    
count = 0
for f in glob.glob("Orthogroup_alignments_trees/trees/*.groupname.treefile"):
    tree = ete3.PhyloTree(f, format=2)
    for l in tree.iter_leaves():
        l.add_feature(pr_name='clade', pr_value=l.name.split('..')[0])
    if find_sisters(tree, 1, 'Rhodelphis'):
        print(f)
        # shutil.copyfile(f.replace('groupname.treefile', 'nex'), "/local/two/Exchange/projects/Picozoa/20_endosymbiotic_origin_genes/orthogroup_EGTs/{}".format(os.path.basename(f).replace('groupname.treefile', 'nex')))
        count += 1
print(count)

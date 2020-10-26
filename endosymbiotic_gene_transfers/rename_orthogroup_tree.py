import pandas as pd
import ete3
import sys

treefile = sys.argv[1] 

focus_taxon = sys.argv[2]

df = pd.read_csv("taxonomy_orthofinder_selection.csv", sep='\t')
tree = ete3.PhyloTree(treefile, format=2)

df.loc[df['taxonomy'].apply(lambda x: focus_taxon in x), 'group'] = focus_taxon

for l in tree.iter_leaves():
    sp = l.name.split('..')[0]
    group = df[df['Name'] == sp].iloc[0]['group']
    group = group.replace(' ', '_').replace('"', '')
    l.name = "{}..{}".format(group, l.name)
    
tree.write(outfile=treefile.replace('.treefile', '.{}.treefile'.format(focus_taxon)), format=2)

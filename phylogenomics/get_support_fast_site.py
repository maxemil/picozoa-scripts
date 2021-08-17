import ete3
import glob
import os 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor
    

def check_monophyly_support(tree, clades):
    prr = [l for l in tree.get_leaves() if any([c in l.name for c in clades])]
    anc = get_ancestor(prr)
    if all([l in prr for l in anc.get_leaves()]):
        return anc.support
    else:
        return 0.0

steps = ["step{}".format(i) for i in range(11)]
clades = ['Picozoa+Rhodelphis+Rhodophyta', 'Rhodelphis+Rhodophyta', 
          'Archaeplastida', "Archaeplastida+Cryptista", 'Amorphea', 'TSAR', 'SAR']
df = pd.DataFrame(index=steps, columns=clades)
        
for f in glob.glob("*.treefile"):
    step = f.replace('.treefile','')
    tree = ete3.PhyloTree(f)
    df.loc[step, 'Picozoa+Rhodelphis+Rhodophyta'] = check_monophyly_support(tree, ['Picozoa', 'Rhodelphis', 'Rhodophyta'])
    df.loc[step, 'Rhodelphis+Rhodophyta'] = check_monophyly_support(tree, ['Rhodelphis', 'Rhodophyta'])
    df.loc[step, 'Archaeplastida'] = check_monophyly_support(tree, ['Picozoa', 'Archaeplastida'])
    df.loc[step, 'Archaeplastida+Cryptista'] = check_monophyly_support(tree, ['Picozoa', 'Archaeplastida', 'Cryptista'])
    df.loc[step, 'Amorphea'] = check_monophyly_support(tree, ['Opisthokonta', 'Amoebozoa', 'Apusomonadida', 'Breviatea'])
    df.loc[step, 'TSAR'] = check_monophyly_support(tree, ['Telonemia', 'Stramenopila', 'Rhizaria', 'Alveolata'])
    df.loc[step, 'SAR'] = check_monophyly_support(tree, ['Stramenopila', 'Rhizaria', 'Alveolata'])


df.index = [(i+1)*5000 for i in range(11)]

fig, ax = plt.subplots(1, 1, figsize=(10,5))
df.plot.line(ax=ax)
ax.set_xlabel('Number of removed sites')
ax.set_ylabel('UFboot support')
ax.set_xticks([(i+1)*5000 for i in range(11)]);
plt.tight_layout()
plt.savefig('fast_site_removal_supports.pdf')
plt.close()

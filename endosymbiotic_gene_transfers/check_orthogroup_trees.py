import glob
import ete3
from collections import defaultdict
import shutil
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import multiprocessing as mp
import copy
from Bio import SeqIO

# import upsetplot
PLASTID_GROUP = ['Apicomplexa',
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
                 'Rhodelphis',
                 'Rhodophyta']


def find_sisters(tree, mincount, taxon, arabidopsis_egts):
    found_sister = False
    plastid_groups = PLASTID_GROUP + [taxon]
    for n in tree.traverse():
        if not n.is_leaf():
            taxa_count_c1, plastid_count_c1 = count_clade(n.get_children()[0], plastid_groups)
            taxa_count_c2, plastid_count_c2 = count_clade(n.get_children()[1], plastid_groups)
            if check_cyano_sister(taxa_count_c1, taxa_count_c2, 
                                  plastid_count_c1, plastid_count_c2, 
                                  taxon, mincount, n.get_children()[0],
                                  plastid_groups):
                # print(get_mono_clades(n, taxon))
                for l in n.get_children()[0].iter_leaves():
                    if l.clade == taxon:
                        arabidopsis_egts.append(l.name)
                found_sister = True
            elif check_cyano_sister(taxa_count_c2, taxa_count_c1, 
                                  plastid_count_c2, plastid_count_c1, 
                                  taxon, mincount, n.get_children()[1],
                                  plastid_groups):
                # print(get_mono_clades(n, taxon))
                for l in n.get_children()[1].iter_leaves():
                    if l.clade == taxon:
                        arabidopsis_egts.append(l.name)
                found_sister = True
    return found_sister

def check_cyano_sister(tc1, tc2, pc1, pc2, taxon, mincount, lca, plastid_groups):
    if taxon == 'Picozoa':
        if not check_Picozoa_mono(lca, taxon):
            return False
    if is_clade_monophyletic(pc1, 'plastid') \
            and tc1[taxon] >= mincount \
            and check_direct_sister(lca, taxon, plastid_groups) \
            and is_clade_monophyletic(tc2, 'Cyanobacteria', impurity=0.15) \
            and tc2['Cyanobacteria'] >= 2:
        return True

def check_Picozoa_mono(lca, taxon):
    if len(set([l.name.split('..')[1] for l in lca.get_leaves() if l.clade == taxon])) < 2:
        return False
    for clade in get_mono_clades(lca, taxon):
        if len(set([l.name.split('..')[1] for l in clade.get_leaves()])) > 1:
            return True
    return False
    
def get_mono_clades(node, taxon):
    seeds = set([l for l in node.get_leaves() if l.clade == taxon])
    nodes = set()
    for s in seeds:
        n = s
        while all([l.clade == taxon for l in n.up.get_leaves()]):
            n = n.up
        nodes.add(n)
    return nodes
                
def check_direct_sister(node, taxon, plastid_groups):
    euk_plastid_groups = copy.copy(plastid_groups)
    if not taxon == 'Paulinella': 
        euk_plastid_groups.remove('Cyanobacteria')
        euk_plastid_groups.remove('Paulinella')
    for n in get_mono_clades(node, taxon):
        if taxon == 'Picozoa':
            if not check_Picozoa_mono(n, taxon):
                continue
        sister = n.up.get_children()[0] if not n.up.get_children()[0] == n else n.up.get_children()[1]
        taxa_count, plastid_count = count_clade(sister, euk_plastid_groups)
        if is_clade_monophyletic(plastid_count, 'plastid', impurity=0.01):
            return True

def count_clade(node, plastid_groups):
    taxa_count = count_taxa(node)
    plastid_count = count_plastid(taxa_count, plastid_groups)
    return taxa_count, plastid_count

def count_plastid(taxa_count, plastid_groups):
    plastid_count = {'plastid':0, 'other':0}
    for k,v in taxa_count.items():
        if k in plastid_groups:
            plastid_count['plastid'] += v
        else:
            plastid_count['other'] += v
    return plastid_count

def is_clade_monophyletic(taxa_count, taxon, impurity=0.1):
    if taxa_count[taxon] >= (1-impurity) * sum(taxa_count.values()):
        return True
    else:
        return False

def count_taxa(node):
    taxa_count = defaultdict(int)
    for l in node.iter_leaves():
        taxa_count[l.clade] += 1
    return taxa_count

def parse_contamination():
    contamination_contigs = []
    for file in glob.glob("../21_contamination_final_assemblies/*contamination.fasta"):
        for rec in SeqIO.parse(file, 'fasta'):
            contamination_contigs.append(rec.id.split('.')[0])
    return contamination_contigs
    
def check_tree(f, taxon):
    arabidopsis_egts = []
    tree = ete3.PhyloTree(f, format=2)
    contamination_contigs = parse_contamination()
    for l in tree.iter_leaves():
        if 'Uroglena_WA34KE' in l.name:
            l.add_feature(pr_name='clade', pr_value='contamination')
        elif l.name.split('..')[2].split('.')[0] in contamination_contigs:
            l.add_feature(pr_name='clade', pr_value='contamination')
        else:
            l.add_feature(pr_name='clade', pr_value=l.name.split('..')[0])
    min_count = 2 if taxon == 'Picozoa' else 1
    if find_sisters(tree, min_count, taxon, arabidopsis_egts):
        # shutil.copyfile(f, "Orthogroup_Selection_{}/EGTs/{}".format(taxon, os.path.basename(f)))
        # shutil.copyfile(f.replace('treefile', 'nex'), "Orthogroup_Selection_{}/EGTs/{}".format(taxon, os.path.basename(f).replace('treefile', 'nex')))
        with open('arabidopsis_egts.csv', 'a') as out:
            for p in arabidopsis_egts:
                print("{}\t{}".format(f, p), file=out)
        return os.path.basename(f).replace('.{}.treefile'.format(taxon), '')
    else:
        tree.set_outgroup(tree.get_midpoint_outgroup())
        if find_sisters(tree, min_count, taxon, arabidopsis_egts):
            # shutil.copyfile(f, "Orthogroup_Selection_{}/EGTs/{}".format(taxon, os.path.basename(f)))
            # shutil.copyfile(f.replace('treefile', 'nex'), "Orthogroup_Selection_{}/EGTs/{}".format(taxon, os.path.basename(f).replace('treefile', 'nex')))
            with open('arabidopsis_egts.csv', 'a') as out:
                for p in arabidopsis_egts:
                    print("{}\t{}".format(f, p), file=out)
            return os.path.basename(f).replace('.{}.treefile'.format(taxon), '')

taxon2count = {}
focus_taxa = ['Rattus',
              'Telonema',
              'Dictyostelium',
              'Tetrahymena',
              'Thecamonas',
              'Phytophthora',
              'Neurospora',
              'Galdieria',
              'Chloropicon',
              'Arabidopsis',
              'Cyanophora',
              'Emiliania',
              'Bigelowiella',
              'Leptocylindrus',
              'Guillardia',
              'Vitrella',
              'Paulinella',
              'Dinobryon_sp_UTEXLB2267',
              'Mallomonas',
              'Ochromonadales_sp_CCMP2298',
              'Alexandrium',
              'Pedospumella_elongata',
              'Paraphysomonas_bandaiensis',
              'Spumella_bureschii_JBL14',
              'Rhodelphis',
              'Toxoplasma',
              'Cryptosporidium',
              'Hematodinium',
              'Helicosporidium',
              'Polytomella', 
              'Goniomonas',
              'Cryptomonas',
              'Picozoa']

pool = mp.Pool(4)
for taxon in focus_taxa:
    os.makedirs("Orthogroup_Selection_{}/EGTs".format(taxon), exist_ok=True)
    files = glob.glob("Orthogroup_Selection_{}/EGTs/*.{}.treefile".format(taxon, taxon))
    results = pool.starmap(check_tree, list(zip(files, [taxon] * len(files))))
    results = [i for i in results if i]
    print(taxon, len(results))
    taxon2count[taxon] = results
pool.close()

from collections import Counter
df = pd.DataFrame({k:Counter(v) for k, v in taxon2count.items()}).T.fillna(0).astype(int)
df.T.to_csv('EGT_OGs.csv', sep='\t', index=True, header=True)

pd.Series(['OG0001421','OG0001259','OG0002552','OG0002725','OG0003800','OG0003816',
            'OG0004004','OG0004116','OG0004388','OG0009417','OG0011992','OG0018782',
            'OG0023785']).apply(lambda x: x in taxon2count['Picozoa']).sum()

tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
taxon2label = {}
for t in focus_taxa:
    taxon2label[t] = tbl.loc[tbl['Name'].apply(lambda x: t in x), 'Name'].iloc[0].replace('_', ' ')
taxon2label['Rhodelphis'] = 'Rhodelphis'
taxon2label['Picozoa'] = 'Picozoa'


cols = {} 
for t in ['Galdieria', 'Chloropicon', 'Arabidopsis', 'Cyanophora', 'Emiliania',
      'Bigelowiella', 'Leptocylindrus', 'Guillardia', 'Vitrella', 'Paulinella',
      'Dinobryon_sp_UTEXLB2267', 'Mallomonas','Ochromonadales_sp_CCMP2298',
      'Alexandrium']:
      cols[t] = 'green' 
for t in ['Pedospumella_elongata', 'Paraphysomonas_bandaiensis',
      'Spumella_bureschii_JBL14', 'Rhodelphis', 'Toxoplasma', 'Helicosporidium',
      'Polytomella', 'Goniomonas', 'Cryptomonas']:
      cols[t] = 'blue'
for t in ['Cryptosporidium', 'Hematodinium']:
    cols[t] = 'darkred'
for t in ['Rattus', 'Telonema', 'Dictyostelium', 'Tetrahymena','Thecamonas', 
       'Phytophthora', 'Neurospora']:
       cols[t] = 'black'
cols['Picozoa'] = 'red'

fig, axs = plt.subplots(1, 2, figsize=(20,10))

df = pd.DataFrame(taxon2count.keys())
df.columns = ['Taxon']
df['Count'] = df['Taxon'].apply(lambda x: len(taxon2count[x]))
df['Color'] = df['Taxon'].apply(lambda x: cols[x])
df.sort_values(by='Count', inplace=True)
df['Label'] = df.apply(lambda x: "{} ({})".format(taxon2label[x['Taxon']], x['Count']), axis=1)
df['Label2'] = df.apply(lambda x: "{}".format(taxon2label[x['Taxon']], x['Count']), axis=1)
sns.barplot(x="Count", y="Label", data=df, palette=df['Color'], ax=axs[0])
axs[0].set_ylabel('Taxon')

for tick in axs[0].get_yticklabels():
    tick.set_color(df.loc[df['Label'] == tick.get_text(), 'Color'].iloc[0])
axs[0].spines["right"].set_visible(False)
axs[0].spines["top"].set_visible(False)
axs[0].spines["bottom"].set_visible(False)
axs[0].spines["left"].set_visible(False)

p2p = pd.DataFrame(index=df['Taxon'], columns=df['Taxon'])
for t1 in df['Taxon']:
    for t2 in df['Taxon']:
        p2p.loc[t1,t2] = len(set(taxon2count[t1]).intersection(set(taxon2count[t2])))
p2p = p2p[p2p.columns].astype(int)
mask = np.zeros_like(p2p)
mask[np.triu_indices_from(mask)] = True
sns.heatmap(p2p, xticklabels=True, yticklabels=True, cmap="YlGnBu", ax=axs[1], mask=mask)
axs[1].set_yticklabels('')
axs[1].set_xticklabels(df['Label2'])
axs[1].set_ylabel('')
for tick in axs[1].get_xticklabels():
    tick.set_color(df.loc[df['Label2'] == tick.get_text(), 'Color'].iloc[0])
plt.tight_layout()
plt.savefig('EGT_groups.pdf')
plt.clf()

# tc = upsetplot.from_contents(taxon2count)
# fig, axs = plt.subplots(figsize=(30,40))
# upsetplot.plot(tc, sort_by='cardinality', fig=fig, show_counts=True)
# plt.savefig('test.pdf')

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

def find_sisters(tree, mincount, taxon, group):
    plastid_groups = [taxon, group]
    for n in tree.traverse():
        if not n.is_leaf():
            taxa_count_c1, plastid_count_c1 = count_clade(n.get_children()[0], plastid_groups)
            taxa_count_c2, plastid_count_c2 = count_clade(n.get_children()[1], plastid_groups)
            if check_cyano_sister(taxa_count_c1, taxa_count_c2, 
                                  plastid_count_c1, plastid_count_c2, 
                                  taxon, mincount, n.get_children()[0],
                                  plastid_groups):
                # print(get_mono_clades(n, taxon))
                return True
            elif check_cyano_sister(taxa_count_c2, taxa_count_c1, 
                                  plastid_count_c2, plastid_count_c1, 
                                  taxon, mincount, n.get_children()[1],
                                  plastid_groups):
                # print(get_mono_clades(n, taxon))
                return True
    return False

def check_cyano_sister(tc1, tc2, pc1, pc2, taxon, mincount, lca, plastid_groups):
    if taxon == 'Picozoa':
        if not check_Picozoa_mono(lca, taxon):
            return False
    if is_clade_monophyletic(pc1, 'plastid') \
            and tc1[taxon] >= mincount \
            and is_clade_monophyletic(tc2, 'Bacteria'):
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

def check_tree(f, taxon, group):
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
    if find_sisters(tree, min_count, taxon, group):
        # shutil.copyfile(f, "Orthogroup_LGT_selections/{}/EGTs/{}".format(taxon, os.path.basename(f)))
        # shutil.copyfile(f.replace('treefile', 'nex'), "Orthogroup_LGT_selections/{}/EGTs/{}".format(taxon, os.path.basename(f).replace('treefile', 'nex')))
        return os.path.basename(f).replace('.{}.treefile'.format(taxon), '')
    else:
        tree.set_outgroup(tree.get_midpoint_outgroup())
        if find_sisters(tree, min_count, taxon, group):
            # shutil.copyfile(f, "Orthogroup_LGT_selections/{}/EGTs/{}".format(taxon, os.path.basename(f)))
            # shutil.copyfile(f.replace('treefile', 'nex'), "Orthogroup_LGT_selections/{}/EGTs/{}".format(taxon, os.path.basename(f).replace('treefile', 'nex')))
            return os.path.basename(f).replace('.{}.treefile'.format(taxon), '')

taxon2count = {}
focus_taxa = [('Picozoa', 'Picozoa'),
              ('Chloropicon', 'Chlorophyta'),
              ('Telonema', 'Telonema'),
              ('Cyanophora', 'Glaucophyta'),
              ('Galdieria', 'Rhodophyta'),
              ('Paulinella', 'Paulinella'),
              ('Rhodelphis', 'Rhodelphis'),
              ('Rattus','Metazoa'),
              ('Dictyostelium','Amoebozoa'),
              ('Tetrahymena','Ciliophora'),
              ('Thecamonas','Apusomonadida'),
              ('Phytophthora','Peronosporomycetes'),
              ('Neurospora','Fungi'),
              ('Arabidopsis','Streptophyta'),
              ('Emiliania','Haptophyta'),
              ('Bigelowiella','Chlorarachniophyceae'),
              ('Leptocylindrus','Diatomeae'),
              ('Guillardia','Cryptophyceae'),
              ('Vitrella','Colpodellida'),
              ('Alexandrium','Dinoflagellata'),
              ('Toxoplasma','Apicomplexa'),
              ('Cryptosporidium','Apicomplexa'),
              ('Hematodinium','Dinoflagellata'),
              ('Helicosporidium','Chlorophyta'),
              ('Goniomonas','Cryptophyceae'),
              ('Cryptomonas','Cryptophyceae'),
              ('Dinobryon_sp_UTEXLB2267','Chrysophyceae'),
              ('Ochromonadales_sp_CCMP2298','Chrysophyceae'),
              ('Pedospumella_elongata','Chrysophyceae'),
              ('Paraphysomonas_bandaiensis','Chrysophyceae'),
              ('Spumella_bureschii_JBL14','Chrysophyceae'),
              ('Mallomonas','Chrysophyceae')]

# focus_taxa = [('Picozoa', 'Picozoa')]


pool = mp.Pool(30)
for taxon, group in focus_taxa:
    os.makedirs("Orthogroup_LGT_selections/{}/EGTs".format(taxon), exist_ok=True)
    files = glob.glob("Orthogroup_LGT_selections/{}/fast_trees/*.{}.treefile".format(taxon, taxon))
    results = pool.starmap(check_tree, list(zip(files, [taxon] * len(files), [group] * len(files))))
    results = [i for i in results if i]
    print(taxon, len(results))
    taxon2count[taxon] = results
pool.close()


pd.Series(['OG0001421','OG0001259','OG0002552','OG0002725','OG0003800','OG0003816',
            'OG0004004','OG0004116','OG0004388','OG0009417','OG0011992','OG0018782',
            'OG0023785']).apply(lambda x: x in taxon2count['Picozoa']).sum()

tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
taxon2label = {}
for t, g in focus_taxa:
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
      'Polytomella', 'Cryptomonas']:
      cols[t] = 'blue'
for t in ['Cryptosporidium', 'Hematodinium']:
    cols[t] = 'darkred'
for t in ['Rattus', 'Telonema', 'Dictyostelium', 'Tetrahymena','Thecamonas', 
       'Phytophthora', 'Neurospora', 'Goniomonas']:
       cols[t] = 'black'
cols['Picozoa'] = 'red'

fig, ax = plt.subplots(1, 1, figsize=(10,5))

df = pd.DataFrame(taxon2count.keys())
df.columns = ['Taxon']
df['Count'] = df['Taxon'].apply(lambda x: len(taxon2count[x]))
df['Color'] = df['Taxon'].apply(lambda x: cols[x])
df.sort_values(by='Count', inplace=True)
df['Label'] = df.apply(lambda x: "{} ({})".format(taxon2label[x['Taxon']], x['Count']), axis=1)
df['Label2'] = df.apply(lambda x: "{}".format(taxon2label[x['Taxon']], x['Count']), axis=1)
sns.barplot(x="Count", y="Label", data=df, palette=df['Color'], ax=ax)
ax.set_ylabel('Taxon')

for tick in ax.get_yticklabels():
    tick.set_color(df.loc[df['Label'] == tick.get_text(), 'Color'].iloc[0])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)


plt.tight_layout()
plt.savefig('LGT_groups.pdf')
plt.close()

df['EGT_count'] = df['Taxon'].apply(lambda x: len(taxon2count[x]))
df['EGT/LGT Ratio'] =  df['EGT_count'] / df['Count']
df.sort_values(by='EGT/LGT Ratio', inplace=True)

fig, ax = plt.subplots(1, 1, figsize=(10,5))
sns.barplot(x="EGT/LGT Ratio", y="Label2", data=df, palette=df['Color'], ax=ax)
ax.set_ylabel('Taxon')

for tick in ax.get_yticklabels():
    tick.set_color(df.loc[df['Label2'] == tick.get_text(), 'Color'].iloc[0])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.tight_layout()
plt.savefig('Ratio_groups.pdf')
plt.close()

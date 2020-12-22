from collections import defaultdict
import glob
import matplotlib.pyplot as plt
import os 
import pandas as pd
import seaborn as sns
from Bio import SeqIO

Phylogenomic = defaultdict(float)
all_Phylogenomic = 0
for f in glob.glob('../14_170_taxa_selection/320_genes/raw_seqs/*'):
    specs = set()
    for rec in SeqIO.parse(f, 'fasta'):
        if rec.id.startswith('Picozoa'):
            spec = rec.id.split('@')[0].split('_')[-1]
            if spec.startswith('SAG') or spec.startswith('COSAG'): 
                specs.add(spec)
    for s in specs:
        Phylogenomic[s] += 1
    if specs:
        all_Phylogenomic += 1
for k, v in Phylogenomic.items():
    Phylogenomic[k] = v/320 * 100
all_Phylogenomic = all_Phylogenomic/320 * 100

BUSCO = defaultdict(float)
for f in glob.glob('../07_BUSCO/eukaryota_odb9_results/run_*/full_table*.tsv'):
    spec = os.path.basename(f).replace("full_table_", "").replace('.tsv', '')
    bo = pd.read_csv(f, sep ='\t', skiprows=4)
    bo = bo.drop_duplicates('# Busco id')
    bo["Status"] = bo["Status"].apply(lambda x: 0 if x== "Missing" else 1)
    assert len(bo['Status']) == 303
    BUSCO[spec] = bo['Status'].sum()/len(bo['Status']) * 100

all_BUSCO = BUSCO['all_assemblies']

df = pd.DataFrame.from_dict(BUSCO, orient='index', columns=['Completeness'])
df['type'] = df.index.map(lambda x: 'SAG' if x.startswith('SAG') else 'COSAG')
df['dataset'] = 'BUSCO'
df = df[df.index.map(lambda x: x in Phylogenomic)]

df2 = pd.DataFrame.from_dict(Phylogenomic, orient='index', columns=['Completeness'])
df2['type'] = df2.index.map(lambda x: 'SAG' if x.startswith('SAG') else 'COSAG')
df2['dataset'] = 'Phylogenomic'
df2 = df2[df2.index.map(lambda x: x in df.index)]
df = df.append(df2)
df = df.append(pd.DataFrame([['combined','BUSCO',all_BUSCO]], columns=['type', 'dataset', 'Completeness'], index=['All']))
df = df.append(pd.DataFrame([['combined','Phylogenomic',all_Phylogenomic]], columns=['type', 'dataset', 'Completeness'], index=['All']))

fig, ax = plt.subplots(1, 1, figsize=(8,5))

sns.boxplot(x="Completeness", y="type", hue='dataset', data=df, ax=ax)
leg = ax.get_legend()
new_labels = ['BUSCO (303)', 'Phylogenomics (320)']
for t, l in zip(leg.texts, new_labels): t.set_text(l)
ax.set(xlim=(0,100))
plt.tight_layout()
plt.savefig('Completeness.pdf')
plt.close()

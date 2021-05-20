import pandas as pd
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

egt_ogs = pd.read_csv('EGT_OGs.csv', sep="\t", index_col=0, header=0)
lgt_ogs = pd.read_csv('LGT_OGs.csv', sep="\t", index_col=0, header=0)


# egt_seqs = defaultdict(lambda:defaultdict(list))
# 
# for line in open('EGT_SEQs.csv'):
#     line = line.strip().split()
#     og = os.path.basename(line[0]).split('.')[0]
#     sp = os.path.basename(line[0]).split('.')[1]
#     egt_seqs[sp][og].append(line[1])
# 
# for sp, og2seqs in egt_seqs.items():
#     assert len(og2seqs.keys()) == egt_ogs[sp].sum()
#     for og, seqs in og2seqs.items():
#         egt_ogs.loc[og,sp] = len(seqs)
# 
# egt_ogs = egt_ogs.fillna(0)




tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
taxon2label = {}
for t in egt_ogs.columns:
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




df = pd.DataFrame(egt_ogs.columns)
df.columns = ['Taxon']
df['EGT'] = df['Taxon'].apply(lambda x: egt_ogs[x].sum())
df['LGT'] = df['Taxon'].apply(lambda x: lgt_ogs[x].sum())
df['Color'] = df['Taxon'].apply(lambda x: cols[x])
df['Label_egt'] = df.apply(lambda x: "{} ({})".format(taxon2label[x['Taxon']], x['EGT']), axis=1)
df['Label_lgt'] = df.apply(lambda x: "{} ({})".format(taxon2label[x['Taxon']], x['LGT']), axis=1)
df['Label'] = df.apply(lambda x: "{}".format(taxon2label[x['Taxon']]), axis=1)
df['EGT/LGT Ratio'] =  df['EGT'] / df['LGT']
# df.to_csv('EGT_LGT_results.csv', sep='\t', header=True, index=False)

fig, ax = plt.subplots(1, 2, figsize=(10,5))
df.sort_values(by='EGT', inplace=True)
sns.barplot(x="EGT", y="Label_egt", data=df, palette=df['Color'], ax=ax[0])
ax[0].set_ylabel('Taxon')
for tick in ax[0].get_yticklabels():
    tick.set_color(df.loc[df['Label_egt'] == tick.get_text(), 'Color'].iloc[0])
ax[0].spines["right"].set_visible(False)
ax[0].spines["top"].set_visible(False)
ax[0].spines["bottom"].set_visible(False)
ax[0].spines["left"].set_visible(False)
# plt.tight_layout()
# plt.savefig('EGT_groups.pdf')
# plt.close()

# fig, ax = plt.subplots(1, 1, figsize=(10,5))
# df.sort_values(by='EGT/LGT Ratio', inplace=True)
sns.barplot(x="EGT/LGT Ratio", y="Label", data=df, palette=df['Color'], ax=ax[1])
ax[1].set_ylabel('Taxon')
for tick in ax[1].get_yticklabels():
    tick.set_color(df.loc[df['Label'] == tick.get_text(), 'Color'].iloc[0])
ax[1].spines["right"].set_visible(False)
ax[1].spines["top"].set_visible(False)
ax[1].spines["bottom"].set_visible(False)
ax[1].spines["left"].set_visible(False)
ax[1].set_xscale("log")
plt.axvline(x=1, color='black', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('EGT_Ratio_groups.pdf')
# plt.savefig('Ratio_groups.pdf')
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(10,5))
df.sort_values(by='LGT', inplace=True)
sns.barplot(x="LGT", y="Label_lgt", data=df, palette=df['Color'], ax=ax)
ax.set_ylabel('Taxon')
for tick in ax.get_yticklabels():
    tick.set_color(df.loc[df['Label_lgt'] == tick.get_text(), 'Color'].iloc[0])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
plt.tight_layout()
plt.savefig('LGT_groups.pdf')
plt.close()

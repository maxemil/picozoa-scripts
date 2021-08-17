from collections import defaultdict
import glob
import matplotlib.pyplot as plt
import os 
import pandas as pd
import seaborn as sns


phylo_assemblies = ['SAG11', 'SAG22', 'SAG25', 'SAG31', 'COSAG01', 'COSAG02', 
                    'COSAG03', 'COSAG04', 'COSAG05', 'COSAG06']
for f in glob.glob('*/run_eukaryota_odb10/full_table.tsv'):
    spec = f.split('/')[0]
    if not spec in phylo_assemblies:
        continue
    bo = pd.read_csv(f, sep ='\t', skiprows=2)
    if not 'detailed_busco' in locals():
        detailed_busco = pd.DataFrame(columns=set(bo['# Busco id']))
    for bid in set(bo['# Busco id']):
        bido = bo[bo['# Busco id'] == bid]
        if bido.shape[0] > 1:
            detailed_busco.loc[spec,bid] = 'Duplicated'
        else:
            detailed_busco.loc[spec,bid] = bido['Status'].item()

combined_busco = pd.Series(index=detailed_busco.columns,dtype=str)
for bid in detailed_busco.columns:
    if 'Complete' in list(detailed_busco[bid]):
        combined_busco.loc[bid] = 'Complete'
    elif 'Duplicated' in list(detailed_busco[bid]):
        combined_busco.loc[bid] = 'Duplicated'
    elif 'Fragmented' in list(detailed_busco[bid]):
        combined_busco.loc[bid] = 'Fragmented'
    else:
        combined_busco.loc[bid] = 'Missing'
combined_busco.value_counts()

count_status = pd.DataFrame(columns=['Missing', 'Duplicated', 'Fragmented', 'Complete'])
for sag in detailed_busco.index:
    count_status = count_status.append(detailed_busco.loc[sag].value_counts())
count_status['index'] = count_status.index
melt_counts = pd.melt(count_status, 
            id_vars='index', 
            value_vars=list(count_status.columns[0:4]), # list of days of the week
            var_name='Status', 
            value_name='Value')
melt_counts = melt_counts.fillna(0)
melt_counts['Dataset'] = 'single'
melt_counts = melt_counts.append(pd.DataFrame({'Value':combined_busco.value_counts().values, 'Status':combined_busco.value_counts().index}))
melt_counts['Value'] = melt_counts['Value'] /255* 100
melt_counts = melt_counts.fillna('All')


fig, ax = plt.subplots(1, 1, figsize=(8,5))

sns.boxplot(x="Status", y="Value", hue='Dataset', data=melt_counts, ax=ax)
leg = ax.get_legend()
new_labels = ['(CO)SAGs', 'Combined']
for t, l in zip(leg.texts, new_labels): t.set_text(l)
ax.set(ylim=(0,100), ylabel='% of BUSCOs', xlabel='Status')
plt.tight_layout()
plt.savefig('Detailed_busco.pdf')
plt.close()

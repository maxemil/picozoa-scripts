import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

df = pd.read_csv('quast_0kb_cosags/report.tsv', sep='\t', index_col=0)
df = df.T
df = df.sort_values(by='Total length (>= 0 bp)', ascending=False)
df.to_csv('cosags_report.csv', sep='\t', header=True, index=True)

dfs = pd.read_csv('quast_0kb_sags/report.tsv', sep='\t', index_col=0)
dfs = dfs.T
dfs = dfs.sort_values(by='Total length (>= 0 bp)', ascending=False)
dfs.to_csv('cosags_report.csv', sep='\t', header=True, index=True)

# df_clean = pd.read_csv('../21_contamination_final_assemblies/quast_results/results_2020_12_07_14_08_55/report.tsv', sep='\t', index_col=0)


fig, ax = plt.subplots(1, 1, figsize=(8,5))
sns.barplot(y=df.index, x=df['Total length (>= 0 bp)'], color='blue')
plt.tight_layout()
plt.savefig('Assembly_sizes.pdf')
plt.close()


df_all = pd.read_csv('../02_renamed_assemblies/quast_all_assemblies/report.tsv', sep='\t', index_col=0)
df_all = df_all.T
df_all = df_all.sort_values(by='Total length (>= 0 bp)', ascending=False)
df_all.index = pd.Series(df_all.index).apply(lambda x: x.split('.')[0])

cosags2sags = defaultdict(list)
for line in open('../scripts/SAGs_COSAGs.tab'):
    line = line.strip().split()
    if line[1].startswith('COSAG'):
        cosags2sags[line[1]].append(line[0])

for k, v in cosags2sags.items():
    print(k)
    # print(df.loc[k, 'Total length (>= 0 bp)'])
    # print(df_all.loc[v, 'Total length (>= 0 bp)'])
    sags = df_all.loc[v, 'Total length (>= 0 bp)'].sum()
    cosag = df.loc[k, 'Total length (>= 0 bp)'].sum()
    print(cosag/sags)

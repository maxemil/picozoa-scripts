import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

df = pd.read_csv('quast_contamination/report.tsv', sep='\t', index_col=0)
df = df.T
df = df.sort_values(by='Total length (>= 0 bp)', ascending=False)

assemblies = ['SAG02', 'SAG03', 'SAG05', 'SAG07', 'SAG10', 'SAG11', 
              'SAG22', 'SAG25', 'SAG31', 'SAG33', 'SAG37', 'COSAG01', 'COSAG02', 
               'COSAG03', 'COSAG04', 'COSAG05', 'COSAG06']

as_length = pd.DataFrame(columns=['length [bp]','SAG','status'])
for ass in assemblies:
    index = "{}.clean".format(ass)
    as_length.loc[index] = [df.loc[index,'Total length (>= 0 bp)'], 
                    ass, "eukaryotic"]
    index = "{}.contamination".format(ass)
    as_length.loc[index] = [df.loc[index,'Total length (>= 0 bp)'], 
                    ass, "prokaryotic/viral"]

as_length['total_length'] = 0
for ass in assemblies:
    total_length = as_length.loc[as_length.index.map(lambda x: True if x.startswith(ass) else False), 'length [bp]'].sum()
    as_length.loc[as_length.index.map(lambda x: True if x.startswith(ass) else False),'total_length'] = total_length

as_length['length [%]'] = as_length['length [bp]']/as_length['total_length'] *100

ax = sns.barplot(x="SAG", y="length [%]", hue="status", data=as_length)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
plt.tight_layout()
plt.savefig('Contamination_final_SAGs.pdf')
plt.close()

import pandas as pd
import os 

sample2sag = {}
for line in open('scripts/rename_samples.tab'):
    line = line.strip().split()
    sample2sag[line[0].split('_L008')[0]] = line[1]

df = pd.read_csv("fastq_stats.csv", header=0, sep='\s+')

df['sag'] = df['file'].apply(lambda x: sample2sag[os.path.basename(x).split('_L008')[0]])
df['num_seqs'] = df['num_seqs'].apply(lambda x: pd.to_numeric(x.replace(',', '')))
df.loc[df['num_seqs'] < 30000, 'sag']

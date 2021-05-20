import ast
from Bio import SeqIO
from collections import defaultdict

def get_taxon_pairs(taxon_pair_file):
    with open(taxon_pair_file, "r") as handle:
        contents = handle.read()
        taxon_pairs = ast.literal_eval(contents)
    return taxon_pairs

seqs = []
for rec in SeqIO.parse('initial_dataset/concatenation/Picozoa_760taxa_317genes_divvier.aln', 'fasta'):
    seqs.append(rec.id.split('@')[0].replace('/', '_'))

merge1 = get_taxon_pairs('164_selection/taxon_pairs_dict.txt')
merge2 = get_taxon_pairs('67_selection/taxon_pairs_dict.txt')
merged_sp = defaultdict(str)
for k,v in merge1.items():
    for sp in k:
        merged_sp[sp.replace('/', '_')] = v.replace('/', '_')
for k,v in merge2.items():
    for sp in k:
        merged_sp[sp.replace('/', '_')] = v.replace('/', '_')
for k,v in merged_sp.items():
    if v in merged_sp:
        merged_sp[k] = merged_sp[v]

selection = [line.strip().replace('/', '_') for line in open('67_selection/taxon_selection.txt')]


with open("taxon_selection.csv", 'w') as out:
    for s in seqs:
        merged = merged_sp[s]
        selected ='No'
        if s in selection:
            selected = 'Yes'
            # selection.remove(s)
        elif merged in selection:
            selected = 'Yes'
            # selection.remove(merged)
        print("{}\t{}\t{}".format(s, merged, selected), file=out)


# check that theres 67 distinct selected genomes/OTUs
import pandas as pd
df = pd.read_csv('taxon_selection.csv', sep='\t', header=None)
df.loc[df[1].isna(), 1] = df.loc[df[1].isna(), 0]
assert len(set(df.loc[df[2] == 'Yes', 1])) == 67

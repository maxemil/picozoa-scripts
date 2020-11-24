import pandas as pd
from Bio import SeqIO

tbl = pd.read_csv('../20_endosymbiotic_origin_genes/taxonomy_orthofinder_selection.csv', sep='\t')
tbl = tbl[tbl['taxonomy'].apply(lambda x: 'Eukaryota' in x)]

recdict = {}
for rec in SeqIO.parse('pr2_version_4.12.0_18S_taxo_long.fasta', 'fasta'):
    recdict[rec.id] = rec

def find_longest_seq(recs, recdict):
    seq = recdict[recs[0]]
    for r in recs:
        if len(recdict[r].seq) > len(seq.seq):
            seq = recdict[r]
    return seq

name2recs = {}
for n in tbl['Name']:
    name2recs[n] = []
    for r in recdict.keys():
        if n in r:
            name2recs[n].append(r)

with open('pr2_all_groups_selection.fasta', 'w') as out:
    for n, recs in name2recs.items():
        if recs:
            seq = find_longest_seq(recs, recdict)
            SeqIO.write(seq, out, 'fasta')
        else:
            print("no recs found for {}".format(n))

# 
# seqdict = {}
# for rec in SeqIO.parse('18S_picozoa_all.aln', 'fasta'):
#     seqdict[rec.id] = rec
# with open('18S_picozoa_all.unique.aln', 'w') as out:
#     for rec in seqdict.values():
#         SeqIO.write(rec, out, 'fasta')

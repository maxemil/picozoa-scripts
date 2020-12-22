from Bio import SeqIO

sample2id = {}
for line in open('../scripts/rename_samples.tab'):
    line = line.strip().split()
    sample2id[line[0].split('-')[-1]] = line[1]

with open('initial_SAGs_18S.fasta', 'w') as out:
    for rec in SeqIO.parse('../03_18S_23S_tree/Picozoa_18S_queries.fasta', 'fasta'):
        for s, i in sample2id.items():
            if rec.id.startswith(s):
                rec.id = rec.id.replace(s, i)
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

from Bio import SeqIO
import pandas as pd 


proteomics = pd.read_csv('Chlamydomonas Ciliary Proteins.csv', sep='\t')
high_confidence_prots = list(proteomics['Protein'])


recs = []
for rec in SeqIO.parse('Creinhardtii_281_v5.6.protein.fa', 'fasta'):
    if rec.id in high_confidence_prots:
        recs.append(rec)
        high_confidence_prots.remove(rec.id)
for rec in SeqIO.parse('Creinhardtii_additional_flagellar_seqs.fasta', 'fasta'):
    recs.append(rec)


with open('Creinhardtii_high_conf_flaggelar.fasta', 'w') as out:
    for rec in recs:
        SeqIO.write(rec, out, 'fasta')

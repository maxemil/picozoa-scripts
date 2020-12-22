from Bio import SeqIO
import glob

for f in glob.glob('*.clean.fasta'):
    assembly = f.split('.')[0]
    clean_contigs = []
    with open(f.replace('fasta', 'faa'), 'w') as out:
        for rec in SeqIO.parse(f, 'fasta'):
            clean_contigs.append(assembly + "_" + rec.id)
        for rec in SeqIO.parse('../12_final_assemblies/{}.faa'.format(assembly), 'fasta'):
            contig = "_".join(rec.id.split('_')[0:7])
            if contig in clean_contigs:
                SeqIO.write(rec, out, 'fasta')

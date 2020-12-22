from Bio import SeqIO
import glob

contamination_contigs = []
for f in glob.glob("../21_contamination_final_assemblies/*.contamination.fasta"):
    for rec in SeqIO.parse(f, 'fasta'):
        contamination_contigs.append(rec.id.split('cov')[0])

for f in glob.glob("initial_dataset/alignments/*"):
    for rec in SeqIO.parse(f, 'fasta'):
        if 'Picozoa' in rec.id:
            if rec.id.split('@')[1].split('cov')[0] in contamination_contigs:
                print(f, rec.id)

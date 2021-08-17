from Bio import SeqIO
import glob
import os

cleanrecs = []
for ass in glob.glob('../21_contamination_final_assemblies/*.clean.fasta'):
    spec = os.path.basename(ass).split('.')[0]
    for rec in SeqIO.parse(ass, 'fasta'):
        contigid = "{}_{}".format(spec, rec.id)
        contigid = "_".join(contigid.split('_')[0:5])
        cleanrecs.append(contigid)

for aln in glob.glob('initial_dataset/alignments/*.aln'):
    for rec in SeqIO.parse(aln, 'fasta'):
        if 'SAG' in rec.id and 'Picozoa' in rec.id:
            seqid = rec.id.replace('Picozoa_Picozoa_N/A_N/A_N/A_N/A_N/A_', '').replace('@', '_')
            seqid = "_".join(seqid.split('_')[0:5])
            if not seqid in cleanrecs:
                print(aln)
                print(rec.id)

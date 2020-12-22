import glob
from Bio import AlignIO

for alnf in glob.glob("*/concat*/*.aln"):
    aln = AlignIO.read(alnf, 'fasta')
    print(alnf, aln.get_alignment_length())

for alnf in glob.glob("*/*.aln"):
    aln = AlignIO.read(alnf, 'fasta')
    print(alnf, aln.get_alignment_length())
    

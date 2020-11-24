from collections import defaultdict
from Bio import SeqIO
import ete3
import argparse

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def is_prokaryote(taxid):
    contaminants = [2, 2157, 10239] # Bacteria, Archaea or Viruses, Taxid of the suppposed con
tamination
    #contaminants = [2759] # Eukaryota
    lineage = ncbi.get_lineage(taxid)
    if any([t in lineage for t in contaminants]):
        return True
    return False

def get_contamination(blast_result):
    seq2contam = defaultdict(list)
    for line in open(blast_result):
        line = line.strip().split('\t')
        if float(line[4]) > 30: # Position of the bitscore of the blast hit
            try:
                hitid = line[2].split(';')[0] 
                seq2contam[line[0]].append(is_prokaryote(hitid))
            except:
                pass
    return seq2contam
    

def filter_contaminant_seqs(seq2contam):
    contam_seqs = []
    unclear_seqs = []
    clean_seqs = []
    for seq, hits in seq2contam.items():
        if all(hits):
            contam_seqs.append(seq)
        elif any(hits):
            unclear_seqs.append(seq)
        else:
            clean_seqs.append(seq)
    return clean_seqs, contam_seqs, unclear_seqs


def clean_fasta_file(fastafile, clean_seqs, contam_seqs, unclear_seqs, outfile):
    with open("{}.fasta".format(outfile), 'w') as out, \
         open("{}.contam.fasta".format(outfile), 'w') as contam, \
         open("{}.unclear.fasta".format(outfile), 'w') as unclear:
        for rec in SeqIO.parse(fastafile, "fasta"):
            if rec.id in clean_seqs:
                SeqIO.write(rec, out, 'fasta')
            elif rec.id in contam_seqs:
                SeqIO.write(rec, contam, 'fasta')
            elif rec.id in unclear_seqs:
                SeqIO.write(rec, unclear, 'fasta')

def main(args):
    seq2contam = get_contamination(args.blast)
    print("Finished parsing blast output")
    clean_seqs, contam_seqs, unclear_seqs = filter_contaminant_seqs(seq2contam)
    print("Finished calculating the contamination coverage")
    clean_fasta_file(args.fasta, clean_seqs, contam_seqs, unclear_seqs, args.output)
    print("Finished writing the fasta files")

parser = argparse.ArgumentParser()
parser.add_argument('--blast', '-b', type=str, help='input blast results file')
parser.add_argument('--fasta', '-f', type=str, help='input fasta file with contigs')
parser.add_argument('--output', '-o', type=str, help='output fasta files basename')
args = parser.parse_args()

if __name__ == '__main__':
    main(args)

from collections import defaultdict
from Bio import SeqIO
import ete3
import argparse
import os
import glob

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def is_prokaryote(taxid, contaminants):
    # contaminants = [10239] # Viruses
    # contaminants = [2759] # Eukaryota
    lineage = ncbi.get_lineage(taxid)
    if any([t in lineage for t in contaminants]):
        return True
    return False

def get_contamination(blast_result, contaminants=[2, 2157, 10239]):
    seq2contam = defaultdict(list)
    for line in open(blast_result):
        line = line.strip().split('\t')
        if float(line[4]) > 30 and line[2]: # Position of the bitscore of the blast hit
            try:
                for hitid in line[2].split(';'):
                    seq2contam[line[0]].append(is_prokaryote(hitid, contaminants))
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

def parse_contigs(file):
    seqdict = {}
    spec = os.path.basename(file).replace('.fasta', '')
    for rec in SeqIO.parse(file, 'fasta'):
        seqdict["{}_{}".format(spec, rec.id.split('.')[0])] = rec
    return seqdict

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

def from_prots_2_contig(clean_seqs, contam_seqs):
    contig2contam = defaultdict(list)
    for seq in clean_seqs:
        contig2contam[seq.split('.')[0]].append(False)
    for seq in contam_seqs:
        contig2contam[seq.split('.')[0]].append(True)
    return contig2contam    

def classify_contigs(contig2contam):
    clean_contigs = []
    contam_contigs = []
    for k, v in contig2contam.items():
        if v.count(True) > 0.6 * len(v):
            contam_contigs.append(k)
        else:
            clean_contigs.append(k)
    return clean_contigs, contam_contigs

def write_contigs(seqdict, contigs, outfile):
    with open(outfile, 'w') as out:
        for c in contigs:
            SeqIO.write(seqdict[c], out, 'fasta')

def write_contigs_inverse(seqdict, contigs, outfile):
    with open(outfile, 'w') as out:
        for v, rec in seqdict.items():
            if not v in contigs:
                SeqIO.write(rec, out, 'fasta')

for file in glob.glob("../12_final_assemblies/*.fasta"):
    seqdict = parse_contigs(file)
    spec = os.path.basename(file).replace('.fasta', '')
    # seq2contam = get_contamination("DIAMOND_NR/{}.faa.out".format(spec), contaminants=[10239])    
    seq2contam = get_contamination("DIAMOND_NR/{}.faa.out".format(spec))    
    clean_seqs, contam_seqs, unclear_seqs = filter_contaminant_seqs(seq2contam)
    contig2contam = from_prots_2_contig(clean_seqs, contam_seqs)
    clean_contigs, contam_contigs = classify_contigs(contig2contam)
    write_contigs(seqdict, contam_contigs, "{}.contamination.fasta".format(spec))
    # write_contigs(seqdict, contam_contigs, "{}.virus.fasta".format(spec))
    write_contigs_inverse(seqdict, contam_contigs, "{}.clean.fasta".format(spec))

# def main(args):
#     seq2contam = get_contamination(args.blast, contamination=[10239])
#     print("Finished parsing blast output")
#     clean_seqs, contam_seqs, unclear_seqs = filter_contaminant_seqs(seq2contam)
#     contig2contam = from_prots_2_contig(clean_seqs, contam_seqs)
#     print("Finished calculating the contamination coverage")
#     clean_fasta_file(args.fasta, clean_seqs, contam_seqs, unclear_seqs, args.output)
#     print("Finished writing the fasta files")
# 
# parser = argparse.ArgumentParser()
# parser.add_argument('--blast', '-b', type=str, help='input blast results file')
# parser.add_argument('--fasta', '-f', type=str, help='input fasta file with contigs')
# parser.add_argument('--output', '-o', type=str, help='output fasta files basename')
# args = parser.parse_args()
# 
# if __name__ == '__main__':
#     main(args)

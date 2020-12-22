from Bio import SeqIO, Seq
from collections import defaultdict
import argparse
import ast
import glob
import sys

def check_compatible(seqs, pos):
    aas = set([seq.seq[pos] for seq in seqs])
    if len(aas) == 1:
        return aas.pop()
    elif len(aas) == 2 and '-' in aas:
        aas.remove('-')
        return aas.pop()
    raise KeyError

def find_longer_seq(seqs):
    longest = len(seqs[0].seq.ungap('-'))
    longest_seq = seqs[0]
    for seq in seqs[1:]:
        if len(seq.seq.ungap('-')) > longest:
            longest = len(seq.seq.ungap('-'))
            longest_seq = seq
    return longest_seq

def select_sequence(seqs):
    try:
        merged_seq = ""
        for i in range(len(seqs[0].seq)):
            merged_seq += check_compatible(seqs, i)
        seqs[0].seq = Seq.Seq(merged_seq)
        for rec in seqs[1:]:
            seqs[0].id += '_{}'.format(rec.id.split('@')[1])
        seqs[0].description = ""
        return seqs[0]
    except KeyError:
        return find_longer_seq(seqs)

def check_taxon_duplicates(fasta, taxon_pairs):
    taxon_list = [t for pair in taxon_pairs.keys() for t in pair]
    taxon_dict = defaultdict(list)
    for rec in SeqIO.parse(fasta, 'fasta'):
        if rec.id.split('@')[0] in taxon_list:
            for key in taxon_pairs.keys():
                for tax in key:
                    if rec.id.split('@')[0] == tax:
                        rec.id = "{}@{}".format(taxon_pairs[key], rec.id.split('@')[1])
                        rec.description = ""
                        taxon_dict[key].append(rec)
        else:
            taxon_dict[rec.id.split('@')[0]].append(rec)
    for k, v in taxon_dict.items():
        if len(v) > 1:
            taxon_dict[k] = [select_sequence(v)]
    return taxon_dict


def get_taxon_pairs(taxon_pair_file):
    with open(taxon_pair_file, "r") as handle:
        contents = handle.read()
        taxon_pairs = ast.literal_eval(contents)
    return taxon_pairs

def print_cleaned_alignment(seqdict, outfile, selection, raw_seqs):
    with open(outfile, 'w') as out:
        for v in seqdict.values():
            if v[0].id.split('@')[0] in selection:
                if not len(v) == 1:
                    print(v)
                else:
                    if v[0].id in raw_seqs:
                        SeqIO.write(raw_seqs[v[0].id], out, 'fasta')
                    else:
                        v[0].seq = v[0].seq.ungap('-')
                        SeqIO.write(v[0], out, 'fasta')
                        

def main(infolder, outfolder, fastafolder, selection, taxon_pairs):
    selection = [l.strip().replace('N_A', 'N/A') for l in open(selection)]
    taxon_pairs = get_taxon_pairs(taxon_pairs)
    for f in glob.glob("{}/*.aln".format(infolder)):
        seqdict = check_taxon_duplicates(f, taxon_pairs)
        raw_seqs = {rec.id:rec for rec in SeqIO.parse(f.replace(infolder, fastafolder).replace(".aln", ".faa"), 'fasta')}
        print_cleaned_alignment(seqdict, f.replace(infolder, outfolder).replace('aln', 'faa'), selection, raw_seqs)


parser = argparse.ArgumentParser()
parser.add_argument('--inalignments', '-i', type=str, help='input folder alignments')
parser.add_argument('--infasta', '-f', type=str, help='input folder raw fasta')
parser.add_argument('--outfolder', '-o', type=str, help='outfolder selections')
parser.add_argument('--selection', '-s', type=str, help='input taxon selection')
parser.add_argument('--pairs', '-p', type=str, help='input taxon pairs to merge')
args = parser.parse_args()


if __name__ == '__main__':
    main(args.inalignments, args.outfolder, args.infasta, args.selection, args.pairs)

from collections import defaultdict
from Bio import SeqIO
import ete3
from Bio import SeqRecord, Seq
import re
import os

ncbi = ete3.ncbi_taxonomy.NCBITaxa()


with open('Rhodelphis_plastid.faa', 'w') as out:
    for sp in ['Rmarinus', 'Rlimneticus']:
        for line in open('{}_plastid.csv'.format(sp)):
            if line.startswith('Process'):
                continue
            line = line.strip().split('\t')
            rec = SeqRecord.SeqRecord(Seq.Seq(line[5]), id='{}_{}'.format(sp, re.sub(r"[\ |/|\.|'|&]", '_', line[1])))
            SeqIO.write(rec, out, 'fasta')
os.system('diamond makedb --in Rhodelphis_plastid.faa --db Rhodelphis_plastid.dmnd  --threads 10')
os.system('diamond blastp -d Rhodelphis_plastid.dmnd -q ../17_plastid_genes/Picozoa_all.faa --quiet --threads 20 -o Picozoa_vs_Rhodelphis_plastid.out --outfmt 6 --ultra-sensitive --top 0.1')

cands = []
for line in open('Picozoa_vs_Rhodelphis_plastid.out'):
    cands.append(line.split()[0])

with open('Picozoa_vs_Rhodelphis_plastid.faa', 'w') as out:
    for rec in SeqIO.parse('../17_plastid_genes/Picozoa_all.faa', 'fasta'):
        if rec.id in cands:
            SeqIO.write(rec, out, 'fasta')


ctps = []
for line in open('Picozoa_vs_Rhodelphis_plastid_summary.targetp2'):
    line = line.split('\t')
    if line[1] == 'cTP':# and float(line[5]) > 0.8:
        ctps.append(line[0])


contaminants=[2, 2157, 10239]
for line in open('Picozoa_vs_Rhodelphis_plastid_NR.out'):
        line = line.strip().split('\t')
        if 'SAG33' in line[0]:
            continue
        if float(line[4]) > 30 and line[2]:
            try:
                for hitid in line[2].split(';'):
                    hitid = int(hitid)
                    lineage = ncbi.get_lineage(hitid)
                    if not any([t in lineage for t in contaminants]) and line[0] in ctps:
                        print(line[0],line[1], ncbi.get_taxid_translator([hitid])[hitid])
            except:
                pass
                

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
    
seq2contam = get_contamination("Picozoa_vs_Rhodelphis_plastid_NR_noTop.out")
clean_seqs, contam_seqs, unclear_seqs = filter_contaminant_seqs(seq2contam)

plast = set()
mito = set()
for line in open('Picozoa_vs_Rhodelphis_plastid_NR_noTop.out'):
    line = line.strip().split('\t')
    if 'chloroplast' in line[-1]:
        plast.add(line[0])
    elif 'plastid' in line[-1]:
        plast.add(line[0])
    elif 'mitochondrial' in line[-1]:
        mito.add(line[0])

abc_trans = set()
tic32 = set()
suf = set()
hem = set()
for line in open('Picozoa_vs_Rhodelphis_plastid.out'):
    line = line.strip().split('\t')
    if 'ABC' in line[1]:
        abc_trans.add(line[0])
    elif 'Tic32' in line[1]:
        tic32.add(line[0])
    elif 'Suf' in line[1]:
        suf.add(line[0])
    elif 'Hem' in line[1]:
        hem.add(line[0])

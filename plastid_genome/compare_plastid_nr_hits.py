from Bio import SeqIO
from collections import defaultdict
import glob
import ete3

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

def get_name(taxid):
    lineage = ncbi.get_lineage(taxid)
    for l in lineage:
        if ncbi.get_rank([l])[l] == 'phylum':
            return ncbi.get_taxid_translator([l])[l] 
    return ncbi.get_taxid_translator([taxid])[taxid]

def is_prokaryote(taxid, contaminants):
    lineage = ncbi.get_lineage(taxid)
    if any([t in lineage for t in contaminants]):
        return True
    return False

def parse_plastid_blast(blastout):
    evalues = defaultdict(lambda: 0)
    pidents = defaultdict(lambda: 0)
    for line in open(blastout):
        line = line.split()
        evalue = float(line[11])
        pident = float(line[2])
        if evalue > evalues[line[0]]:
            evalues[line[0]] = evalue
        if pident > pidents[line[0]]:
            pidents[line[0]] = pident
    return (evalues, pidents)
            
def parse_nr_blast(blastout, filter):
    evalues = defaultdict(lambda: 0)
    pidents = defaultdict(lambda: 0)
    tax = defaultdict(lambda: set())
    for line in open(blastout):
        line = line.split('\t')
        if line[0] in filter:
            evalue = float(line[12])
            pident = float(line[4])
            if evalue > evalues[line[0]]:
                evalues[line[0]] = evalue
            if pident > pidents[line[0]]:
                pidents[line[0]] = pident
            try:
                for hitid in line[2].split(';'):
                    tax[line[0]].add(get_name(int(hitid)))
            except:
                pass
    return (evalues, pidents, tax)

for f in glob.glob("*SAG*.out"):
    evalues, pidents = parse_plastid_blast(f)
    nr_blast = "../../21_contamination_final_assemblies/DIAMOND_NR/{}.faa.out".format(f.replace(".out", ''))
    evalues_nr, pidents_nr, tax_nr = parse_nr_blast(nr_blast, list(evalues.keys()))
    for k, v in evalues.items():
        if v > evalues_nr[k] and pidents[k] > pidents_nr[k]:
            print(k, v, pidents[k], evalues_nr[k], pidents_nr[k])

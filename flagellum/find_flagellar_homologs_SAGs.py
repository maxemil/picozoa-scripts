import glob
from Bio import SeqIO


# for faa in ../27_assemblies_share/final_SAGs_COSAGs/*.faa; 
# do 
#   diamond blastp -d Creinhardtii_high_conf_flaggelar.fasta -q $faa \
#         --quiet --threads 20 -o $(basename ${faa%%.faa})_flagellum.out \
#         --outfmt 6 --ultra-sensitive --evalue 1e-15; 
# done

ref_prots = [rec.id for rec in SeqIO.parse('Creinhardtii_high_conf_flaggelar.fasta', 'fasta')]
ref_prots = set(ref_prots)

hit_refs = set()
for f in glob.glob("*_flagellum.out"):
    for line in open(f):
        line = line.split('\t')
        hit_refs.add(line[1])

hit_refs = set()
for f in glob.glob("*_flagellum_pazour.out"):
    for line in open(f):
        line = line.split('\t')
        hit_refs.add(line[1])

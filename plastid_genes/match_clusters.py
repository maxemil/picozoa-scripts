from Bio import SeqIO
import glob
import os
import shutil 

rec2genome = {}
for rec in SeqIO.parse('Rhodelphis_plastid.faa', 'fasta'):
    rec2genome[rec.id] = []
    

for line in open('Rlimneticus_vs_Rhodelphis_plastid.out'):
    line = line.strip('\n').split('\t')
    if float(line[2]) > 99:
        rec2genome[line[0]].append(line[1])
for line in open('Rmarinus_vs_Rhodelphis_plastid.out'):
    line = line.strip('\n').split('\t')
    if float(line[2]) > 99:
        rec2genome[line[0]].append(line[1])

Rhodelphis_plastid = [v[0] for v in rec2genome.values() if v] 

ep2annot = {}
for k,v in rec2genome.items():
    for ep in v:
        ep2annot[ep] = k

count = 0
for k,v in rec2genome.items():
    if len(v) != 1:
        pass
    else:
        count += 1

clsts = set()
for f in glob.glob('../20_endosymbiotic_origin_genes/Orthogroup_Sequences/*'):
    for rec in SeqIO.parse(f, 'fasta'):
        if rec.id in Rhodelphis_plastid:
            clsts.add(os.path.basename(f).replace('.fa', ''))

for c in clsts:
    if len(glob.glob('../20_endosymbiotic_origin_genes/Orthogroup_LGT_trees/trees/{}.treefile'.format(c))) > 0:
        shutil.copyfile('../20_endosymbiotic_origin_genes/Orthogroup_LGT_trees/trees/{}.treefile'.format(c), 'Orthogroup_trees/{}.treefile'.format(c))
    elif len(glob.glob('../20_endosymbiotic_origin_genes/Orthogroup_alignments_trees/trees/{}.treefile'.format(c))) > 0:
        shutil.copyfile('../20_endosymbiotic_origin_genes/Orthogroup_alignments_trees/trees/{}.treefile'.format(c), 'Orthogroup_trees/{}.treefile'.format(c))

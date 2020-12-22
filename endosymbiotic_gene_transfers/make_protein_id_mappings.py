import glob
import os
import pandas as pd
from Bio import SeqIO

for f in glob.glob('EukProt/Chrysophyte_proteins/*.fasta'):
    species = f.replace('EukProt/Chrysophyte_proteins/', '').replace('.fasta', '')
    orthofinder_file = "fasta_orthofinder/{}.fasta".format(species)
    mapping_file = "protein_maps_orthofinder/{}.map".format(species)
    with open(mapping_file, 'w') as out:
        for rec, rec2 in zip(SeqIO.parse(f, 'fasta'), SeqIO.parse(orthofinder_file, 'fasta')):
            assert rec.seq == rec2.seq
            mapid = rec.id.replace('@', '')
            mapid = "{}@{}".format(species, mapid)
            print("{}\t{}".format(mapid, rec.description), file=out)
 
for f in glob.glob('PhyloSkeleton/EGT/genomicData/*/*.faa'):
    species = os.path.dirname(f.replace('PhyloSkeleton/EGT/genomicData/', ''))
    orthofinder_file = "fasta_orthofinder/{}.fasta".format(species)
    mapping_file = "protein_maps_orthofinder/{}.map".format(species)
    with open(mapping_file, 'w') as out:
        for rec, rec2 in zip(SeqIO.parse(f, 'fasta'), SeqIO.parse(orthofinder_file, 'fasta')):
            assert rec.seq == rec2.seq
            mapid = "{}@{}".format(species, rec.id)
            print("{}\t{}".format(mapid, rec.description), file=out)

df = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
df = df[df['ID'].apply(lambda x: x.startswith('EP00'))]
for row in df.iterrows():
    ep = row[1]['ID']
    sp = row[1]['Name']
    f = "EukProt/proteins/{}_{}.fasta".format(ep, sp)
    orthofinder_file = "fasta_orthofinder/{}.fasta".format(sp)
    mapping_file = "protein_maps_orthofinder/{}.map".format(sp)
    count = 0
    with open(mapping_file, 'w') as out:
        for rec, rec2 in zip(SeqIO.parse(f, 'fasta'), SeqIO.parse(orthofinder_file, 'fasta')):
            assert rec.seq == rec2.seq
            mapid = "{}@{}_{}".format(sp, ep, count)
            count += 1
            print("{}\t{}".format(mapid, rec.description), file=out)

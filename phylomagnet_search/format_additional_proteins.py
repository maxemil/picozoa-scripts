from Bio import SeqIO
import pandas as pd
import ete3

ncbi = ete3.ncbi_taxonomy.NCBITaxa()
df = pd.read_csv('../20_endosymbiotic_origin_genes/taxonomy_orthofinder_selection.csv', sep='\t')
df = df[df['taxonomy'].apply(lambda x: 'RR clade' in x)]

with open('Rhodelphis_Rhodophytes.faa', 'w') as out:
    for id, sp in zip(df.ID, df.Name):
        file = '../20_endosymbiotic_origin_genes/EukProt/proteins/{}_{}.fasta'.format(id, sp)
        etesp = sp.replace('_', ' ') if len(sp.split('_')) == 2 else sp.split('_')[0]
        taxid = ncbi.get_name_translator([etesp])[etesp][0]
        for rec in SeqIO.parse(file, 'fasta'):
            rec.id = rec.id.replace('|', '_').replace('.', '_')
            rec.id = "{}.{}".format(taxid, rec.id)
            rec.description = ""
            SeqIO.write(rec, out, 'fasta') 
        

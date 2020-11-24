import pandas as pd
import glob
from Bio import SeqIO
import os 
import shutil
from collections import defaultdict

tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')

focus_taxa = ['Rhodelphis',
              'Picozoa',
              'Galdieria',
              'Chloropicon',
              'Cyanophora',
              'Paulinella',
              'Telonema']
              
for focus_taxon in focus_taxa:
    tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
    tbl.loc[tbl['taxonomy'].apply(lambda x: focus_taxon in x), 'group'] = focus_taxon
    tblgroups = set(tbl['group'])
    
    selection = []
    for f in glob.glob('Orthogroup_Sequences/*'):
        groups = defaultdict(int)
        for rec in SeqIO.parse(f, 'fasta'):
            sp = rec.id.split('@')[0]
            group = tbl[tbl['Name'] == sp].iloc[0]['group']
            groups[group] += 1
        if groups['Bacteria'] > 1 and groups[focus_taxon] > 1 \
            and sum(groups.values()) < 1000:
           selection.append(f)
    os.makedirs('Orthogroup_LGT_selections/{}/'.format(focus_taxon))
    for f in selection:
        shutil.copyfile(f, 'Orthogroup_LGT_selections/{}/{}'.format(focus_taxon, os.path.basename(f)))

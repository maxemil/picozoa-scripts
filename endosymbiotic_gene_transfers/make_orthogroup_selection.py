import pandas as pd
import glob
from Bio import SeqIO
import os 
import shutil
from collections import defaultdict

tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')

cyano_ogs = []
for f in glob.glob('Orthogroup_Sequences/*'):
    groups = defaultdict(int)
    for rec in SeqIO.parse(f, 'fasta'):
        sp = rec.id.split('@')[0]
        group = tbl[tbl['Name'] == sp].iloc[0]['group']
        groups[group] += 1
    if groups['Cyanobacteria'] > 0:
        cyano_ogs.append(f)

focus_taxa = ['Rhodelphis',
                'Amoebozoa',
                'Hematodinium',
                'Cryptosporidium',
                'Polytomella',
                'Helicosporidium',
                'Glaucophyta',
                'Pedospumella',
                'Poteriospumella',
                'Paraphysomonas',
                'Cornospumella',
                'Uroglena',
                'Dinobryon', 
                'Poterioochromonas', 
                'Synura',
                'Telonema',
                'Paulinella']

for focus_taxon in focus_taxa:
    tbl = pd.read_csv('taxonomy_orthofinder_selection.csv', sep='\t')
    tbl.loc[tbl['taxonomy'].apply(lambda x: focus_taxon in x), 'group'] = focus_taxon
    tblgroups = set(tbl['group'])

    selection = []
    for f in cyano_ogs:
        groups = {g:0 for g in tblgroups}
        for rec in SeqIO.parse(f, 'fasta'):
            sp = rec.id.split('@')[0]
            group = tbl[tbl['Name'] == sp].iloc[0]['group']
            groups[group] += 1
        if groups[focus_taxon] > 0 and groups['Cyanobacteria'] > 0 \
            and (groups['Chloroplastida'] > 0 or groups['Rhodophyta'] > 0)  \
            and sum(groups.values()) < 1000:
            selection.append(f)

    os.makedirs('Orthogroup_Selection_{}/'.format(focus_taxon))
    for f in selection:
        shutil.copyfile(f, 'Orthogroup_Selection_{}/{}'.format(focus_taxon, os.path.basename(f)))


os.makedirs('Orthogroup_added_Selection/')
added_ogs = set()
picozoa_ogs = set([os.path.basename(p) for p in glob.glob('Orthogroup_Selection_Picozoa/*')])
for focus_taxon in focus_taxa:
    taxon_ogs = set([os.path.basename(p) for p in glob.glob('Orthogroup_Selection_{}/*'.format(focus_taxon))])
    print('Picozoa ({}) - {} ({}): {}'.format(len(picozoa_ogs), focus_taxon, len(taxon_ogs), len(picozoa_ogs.intersection(taxon_ogs))))
    added_ogs = added_ogs.union(taxon_ogs.difference(picozoa_ogs))
for f in added_ogs:
        # print('Orthogroup_Sequences/{}'.format(f))
        shutil.copyfile('Orthogroup_Sequences/{}'.format(f), 'Orthogroup_added_Selection/{}'.format(f))

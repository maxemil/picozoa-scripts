import pandas as pd
import ete3
import glob
from collections import defaultdict
from Bio import SeqIO
import os 

groups = ['Cyanobacteria',
          'Peronosporomycetes',
          'Amoebozoa',
          'Ancoracysta',
          'Ancyromonadida',
          'Apicomplexa',
          'Apusomonadida',
          'Sagenista',
          'Breviatea',
          'Centroplasthelida',
          'Chloroplastida',
          'Choanoflagellata',
          'Ciliophora',
          'Collodictyonidae',
          'Colpodellida',
          'Colponemidae',
          'Cryptophyceae',
          'Dinoflagellata',
          'Diplonemea',
          'Ochrophyta',
          'Euglenida',
          'Filasterea',
          'Fungi',
          'Glaucophyta',
          'Haptophyta',
          'Hemimastigophora',
          'Heterolobosea',
          'Ichthyosporea',
          'Kathablepharidacea',
          'Malawimonadidae',
          'Mantamonas',
          'Metamonada',
          'Metazoa',
          'Palmophyllophyceae',
          'Palpitomonas',
          'Perkinsea',
          'Pluriformea',
          'Chlorarachniophyceae',
          'Paulinella',
          'Endomyxa',
          'Rhodelphis',
          'Rhodophyta',
          'Rigifilida',
          'Rotosphaerida',
          'Telonemia',
          'Tsukubamonadida',
          'Jakobida',
          'Picozoa']

tax_chrysopyceae = {'Acrispumella_msimbazaensis_JBAF33':"Ochromonadales;Acrispumella;",
 'Cornospumella_fuchlensis_AR4D6':"Ochromonadales;Cornospumella;",
 'Dinobryon_FU22KAK':"Ochromonadales;Dinobryon;",
 'Dinobryon_LO226KS':"Ochromonadales;Dinobryon;",
 'Epipyxis_PR26KG':"Ochromonadales;Dinobryon;Epipyxis;",
 'Ochromonas_LO244K-D':"Ochromonadales;Ochromonas;",
 'Pedospumella_sinomuralis_JBCS23':"Ochromonadales;Pedospumella;Pedospumella sinomuralis;",
 'Poterioochromonas_DS':"Ochromonadales;Poterioochromonas;",
 'Poteriospumella_lacustris_JBNZ41':"Ochromonadales;Poteriospumella;",
 'Spumella_NIES48':"Ochromonadales;",
 'Spumella_bureschii_JBL14':'Ochromonadales;"core Spumella";',
 'Spumella_elongata_CCAP_9551':"Ochromonadales;Pedospumella;Pedospumella elongata;",
 'Spumella_lacusvadoi_JBNZ3':'Ochromonadales;"core Spumella";',
 'Spumella_vulgaris':'Ochromonadales;"core Spumella";',
 'Synura_LO234KE':"Synurales;Synura;",
 'Uroglena_WA34KE':"Ochromonadales;Uroglena;"}


def get_group(taxostring, domain):
    for g in groups:
        if g in taxostring:
            return g
    if domain == 'Bacteria':
        return domain
    else:
        print(taxostring)
    
taxa = set()
for f in glob.glob("fasta_orthofinder/*") + glob.glob("fasta_Picozoa/*"):
    for rec in SeqIO.parse(f, 'fasta'):
        taxa.add(rec.id.split('@')[0])
        break

df = pd.read_csv('EukProt/EukProt_included_data_sets.v02.2020_06_30.txt', sep='\t')
ps = pd.read_csv('PhyloSkeleton/EGT/bacteria.selection.tab', sep='\t')

BUSCO = defaultdict(float)
for d in glob.glob('BUSCO_results_eukaryota/EP*'):
    f = "{}/run_eukaryota_odb10/full_table.tsv".format(d)
    bo = pd.read_csv(f, sep ='\t', skiprows=2)
    bo = bo.drop_duplicates('# Busco id')
    bo["Status"] = bo["Status"].apply(lambda x: 0 if x== "Missing" else 1)
    assert len(bo['Status']) == 255
    BUSCO[d.split('/')[1].split('_')[0]] = bo['Status'].sum()/len(bo['Status'])
df['Completeness'] = df['EukProt_ID'].apply(lambda x: BUSCO[x])

tbl = pd.DataFrame(columns=['Name', 'ID', 'group', 'domain', 'completeness', 'taxonomy'])

for t in taxa:
    euk = df[df['Name_to_Use'] == t]
    bac = ps[ps['shortName'].apply(lambda x: t == x)]
    if len(euk) == 1:
        group = get_group(euk.iloc[0]['Taxonomy_UniEuk'], 'Eukaryota')
        tbl = tbl.append({'Name': euk.iloc[0]['Name_to_Use'],
            'ID': euk.iloc[0]['EukProt_ID'],
            'taxonomy': euk.iloc[0]['Taxonomy_UniEuk'], 
            'group': group, 
            'domain': 'Eukaryota', 
            'completeness': euk.iloc[0]['Completeness']},
            ignore_index=True)
    elif len(bac) == 1:
        group = get_group(bac.iloc[0]['taxostring'], 'Bacteria')
        tbl = tbl.append({'Name': bac.iloc[0]['shortName'],
            'ID': bac.iloc[0]['assembly'],
            'taxonomy': bac.iloc[0]['taxostring'], 
            'group': group, 
            'domain': 'Bacteria', 
            'completeness': bac.iloc[0]['completeness']},
            ignore_index=True)
    elif 'SAG' in t:
        tbl = tbl.append({'Name': t,
            'ID': t,
            'taxonomy': 'Eukaryota;Diaphoretickes;Picozoa;{}'.format(t), 
            'group': 'Picozoa', 
            'domain': 'Eukaryota', 
            'completeness': '-'},
            ignore_index=True)
    elif t in tax_chrysopyceae:
        tbl = tbl.append({'Name': t,
            'ID': t,
            'taxonomy': 'Eukaryota;Diaphoretickes;Sar;Stramenopiles;Gyrista;Ochrophyta;"CS clade";Chrysophyceae;{}{}'.format(tax_chrysopyceae[t], t), 
            'group': 'Ochrophyta', 
            'domain': 'Eukaryota', 
            'completeness': '-'},
            ignore_index=True)
    else:
        print(t)

for d in glob.glob('BUSCO_results_eukaryota/[!busco]*'):
    if not os.path.basename(d).startswith('EP'):
        f = "{}/run_eukaryota_odb10/full_table.tsv".format(d)
        bo = pd.read_csv(f, sep ='\t', skiprows=2)
        bo = bo.drop_duplicates('# Busco id')
        bo["Status"] = bo["Status"].apply(lambda x: 0 if x== "Missing" else 1)
        assert len(bo['Status']) == 255
        tbl.loc[tbl['Name'] == d.split('/')[1], 'completeness'] = bo['Status'].sum()/len(bo['Status'])

tbl.to_csv('taxonomy_orthofinder_selection.csv', sep='\t', index=False, header=True)

picozoa_taxa = set(tbl[tbl['group'] == 'Picozoa']['Name'])
cyano_taxa = set(tbl[tbl['group'] == 'Cyanobacteria']['Name'])

tblgroups = set(tbl['group'])

og2groups = defaultdict(lambda: {g:0 for g in tblgroups})

count = 0
sizes = {}
for f in glob.glob('Orthogroup_Sequences/*'):
    seq_count = 0
    for rec in SeqIO.parse(f, 'fasta'):
        seq_count += 1
        sp = rec.id.split('@')[0]
        group = tbl[tbl['Name'] == sp].iloc[0]['group']
        og2groups[f][group] += 1
    if og2groups[f]['Cyanobacteria'] and og2groups[f]['Picozoa']:
        count += 1
        sizes[f] = seq_count

ogs = pd.DataFrame.from_dict(og2groups, orient='index')

selection = ogs[(ogs['Picozoa'] > 0) & \
                (ogs['Cyanobacteria'] > 0) & \
                ((ogs['Chloroplastida'] > 0) | (ogs['Rhodophyta'] > 0)) & \
                (ogs.sum(axis=1) < 1000)]
import shutil
import os

os.makedirs('Orthogroup_Selection/')
for f in selection.index:
    shutil.copyfile(f, 'Orthogroup_Selection/{}'.format(os.path.basename(f)))

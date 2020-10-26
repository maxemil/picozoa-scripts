import pandas as pd
import ete3
import glob
from collections import defaultdict
from Bio import SeqIO
import os 

def subselect(df, n=2):
    df = df.sort_values('Completeness', ascending=False)
    sel = pd.DataFrame(columns=df.columns)
    for row in df.iterrows():
        if len(sel) >= n:
            break
        if not row[1]['Genus_UniEuk'] in set(sel['Genus_UniEuk']):
            sel = sel.append(row[1], ignore_index=True)
    return sel

df = pd.read_csv('EukProt/EukProt_included_data_sets.v02.2020_06_30.txt', sep='\t')

photo_groups = {
                # 'Chloroplastida':'genus', 
                '"core raphidophytes"':2,
                '"PCR clade"':2,
                '"PS clade"':2,
                'Apicomplexa':'genus', 
                'Bolidophyceae':2,
                'Chlorarachniophyceae':'genus',
                'Chlorophyta':'genus', 
                'Chrysophyceae':'all',
                'Colpodellida':'all', 
                'Colponemidae':'all',
                'Cryptophyceae':'genus', 
                'Diatomeae':10,
                'Dictyochales':2,
                'Dinoflagellata':10, 
                'Euglenida':'all', 
                'Eustigmatales':2,
                'Florenciellales':2,
                'Glaucophyta':'all', 
                'Haptophyta':'genus', 
                'Palmophyllophyceae':'all', 
                'Pelagomonadales':2,
                'Phaeomonas':2,
                'Phaeophyceae':2,
                'Pinguiococcus':2,
                'Rhodelphis':'all', 
                'Rhodophyta':'genus',
                'Sarcinochrysidales':2,
                'Streptophyta':'genus',
                'Synchromophyceae':2,
                'Xanthophyceae':2,
                'Perkinsea':'all'}

BUSCO = defaultdict(float)
for d in glob.glob('BUSCO_results_eukaryota/EP*'):
    f = "{}/run_eukaryota_odb10/full_table.tsv".format(d)
    bo = pd.read_csv(f, sep ='\t', skiprows=2)
    bo = bo.drop_duplicates('# Busco id')
    bo["Status"] = bo["Status"].apply(lambda x: 0 if x== "Missing" else 1)
    assert len(bo['Status']) == 255
    BUSCO[d.split('/')[1].split('_')[0]] = bo['Status'].sum()/len(bo['Status'])

df['Completeness'] = df['EukProt_ID'].apply(lambda x: BUSCO[x])


selection = pd.DataFrame(columns=df.columns)
for sg in set(df['Supergroup_UniEuk']):
    sgf = df[df['Taxonomy_UniEuk'].apply(lambda x: sg in x and not 
                                        any(p in x for p in photo_groups))]
    if len(sgf) > 2:
        selection = selection.append(subselect(sgf, n=2))
    else:
        selection = selection.append(sgf)

for p, num in photo_groups.items():
    pf = df[df['Taxonomy_UniEuk'].apply(lambda x: p in x)]
    if num == 'all':
        selection = selection.append(pf)
    elif num == 'genus':
        selection = selection.append(subselect(pf, n=len(pf)))
    else:
        selection = selection.append(subselect(pf, n=num)) 
                                


tree = ete3.Tree('Eukaryota;')
for tax in df['Taxonomy_UniEuk']:
    parent = tree & "Eukaryota"
    for t in tax.split(';'):
        try:
            parent = tree & t
        except:
            parent.add_child(name=t)
            parent = tree & t            
tree.ladderize()

ts = ete3.TreeStyle()
ts.show_leaf_name = False

def my_layout(node):
    if node.is_leaf() and any([node.name == x.split(';')[-1] for x in selection['Taxonomy_UniEuk']]):
        F = ete3.TextFace(node.name, tight_text=True, fgcolor='red')
    else:
        F = ete3.TextFace(node.name, tight_text=True)
    ete3.add_face_to_node(F, node, column=0, position="branch-right")

ts.layout_fn = my_layout
tree.render('test.pdf', tree_style=ts)

# for f in glob.glob('EukProt/Chrysophyte_proteins/*.fasta'):
#     species = f.replace('EukProt/Chrysophyte_proteins/', '').replace('.fasta', '')
#     outfile = "fasta_orthofinder/{}.fasta".format(species)
#     with open(outfile, 'w') as out:
#         for rec in SeqIO.parse(f, 'fasta'):
#             rec.id = rec.id.replace('@', '')
#             rec.id = "{}@{}".format(species, rec.id)
#             rec.description = ""
#             SeqIO.write(rec, out, 'fasta')
# # 
# for f in glob.glob('PhyloSkeleton/EGT/genomicData/*/*.faa'):
#     species = os.path.dirname(f.replace('PhyloSkeleton/EGT/genomicData/', ''))
#     outfile = "fasta_orthofinder/{}.fasta".format(species)
#     with open(outfile, 'w') as out:
#         for rec in SeqIO.parse(f, 'fasta'):
#             rec.id = "{}@{}".format(species, rec.id)
#             rec.description = ""
#             SeqIO.write(rec, out, 'fasta')
# 
# for row in selection.iterrows():
#     ep = row[1]['EukProt_ID']
#     sp = row[1]['Name_to_Use']
#     with open("fasta_orthofinder/{}.fasta".format(sp), 'w') as out:
#         count = 0
#         for rec in SeqIO.parse("EukProt/proteins/{}_{}.fasta".format(ep, sp), 'fasta'):
#             rec.id = "{}@{}_{}".format(sp, ep, count)
#             count += 1
#             rec.description = ""
#             SeqIO.write(rec, out, 'fasta')

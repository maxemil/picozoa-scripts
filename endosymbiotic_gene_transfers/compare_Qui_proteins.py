import glob
import ete3
import os
from Bio import Entrez

def get_acc_for_gi():
    Entrez.email = 'max-emil.schon@icm.uu.se'
    qui_gi = [l.strip() for l in open('gi_list.csv')]
    gi2acc = {}
    for i in range(0, len(qui_gi), 200):
        gi_fetch = ','.join(qui_gi[i:i+200])
        handle = Entrez.efetch(db="protein", id=gi_fetch, rettype="acc")
        acc_list = handle.read().strip().split()
        for gi, acc in zip(qui_gi[i:i+200], acc_list):
            gi2acc[gi] = acc
    return gi2acc

def parse_protein_map(infile):
    p2id = {}
    for l in open(infile):
        l = l.split()
        p2id[l[1]] = l[0]
    return p2id

def get_overlap(gi2acc, protmap):
    qui_egt_prots = []
    for acc in gi2acc.values():
        if acc in protmap:
            qui_egt_prots.append(protmap[acc])
        else:
            print(acc)
    return qui_egt_prots

def format_prot_ids(protid, prefix):
    return "{}..{}".format(prefix, protid.replace('@', '..'))

def get_predicted_egt(infolder, prefix):
    predicted_egt = []
    for tree in glob.glob('{}/*.treefile'.format(infolder)):
        tree = ete3.PhyloTree(tree) 
        for l in tree.iter_leaves():
            if l.name.startswith(prefix):
                predicted_egt.append(l.name)
    return predicted_egt

def get_accessions_from_prots(prots, acc2prot, prefix):
    prot2acc = {format_prot_ids(v, prefix):k for k,v in acc2prot.items()}
    return [prot2acc[p] for p in prots]

def get_subcellular_location(acc):
    req = urllib.request.Request("https://www.uniprot.org/uniprot/{}.xml".format(acc))
    with urllib.request.urlopen(req) as f:
        rec = SeqIO.read(f, 'uniprot-xml')
        if 'comment_subcellularlocation_location' in rec.annotations:
            loc = rec.annotations['comment_subcellularlocation_location'][0]
            if any([p in loc for p in plastid_key]):
                return 'Plastid'
            elif any([p in loc for p in mito_key]):
                return 'Mitochondrial'

def find_subcellular_location(ids):
    name2loc = defaultdict(str)
    for acc in ids:
        try:
            req = urllib.request.Request("https://www.uniprot.org/uniprot/?query={}&format=tab".format(acc))
            with urllib.request.urlopen(req) as f:
                res = f.read().decode('utf-8').strip()
                if res:
                    res = [x.split('\t') for x in res.split('\n')]
                    df = pd.DataFrame(res[1:], columns=res[0])
                    # if df['Status'][0] == 'reviewed':
                    loc = get_subcellular_location(df['Entry name'][0])
                    if loc:
                        name2loc[l] = loc                        
        except Exception as e: 
            pass
    return name2loc

def get_records(accessions):
    Entrez.email = 'max-emil.schon@icm.uu.se'
    records = []
    rec2locus = {}
    rec2organelle = {}
    for i in range(0, len(accessions), 200):
        acc = ','.join(accessions[i:i+200])
        handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode='text')
        # acc_list = handle.read().strip().split()
        for rec in SeqIO.parse(handle, 'genbank'):
            records.append(rec)
            for f in rec.features:
                if f.type == 'CDS':
                    rec2locus[rec.id] = f.qualifiers['locus_tag'][0]
                if f.type == 'source':
                    if 'organelle' in f.qualifiers:
                        rec2organelle[rec.id] = f.qualifiers['organelle'][0]
    return (records, rec2locus, rec2organelle)

egt_trees = [os.path.basename(f).split('.')[0] for f in glob.glob('../Orthogroup_Selection_Chloropicon/EGTs/*.treefile')]

unrecognized_egt = []
for treefile in glob.glob('../Orthogroup_Selection_Chloropicon/trees/*Chloropicon.treefile'):
    if not os.path.basename(treefile).split('.')[0] in egt_trees:
        tree = ete3.PhyloTree(treefile)
        for l in tree.iter_leaves():
            if any([l.name.split('..')[2] in p for p in qui_egt_crei]):
                unrecognized_egt.append(treefile)


gi2acc = get_acc_for_gi()

a_thaliana = parse_protein_map('../protein_maps_orthofinder/Arabidopsis_thaliana.map')
qui_egt_prots = get_overlap(gi2acc, a_thaliana)
qui_egt_prots = [format_prot_ids(p, 'Arabidopsis') for p in qui_egt_prots]
max_egt_prots = get_predicted_egt('../Orthogroup_Selection_Arabidopsis/EGTs', 'Arabidopsis..Arabidopsis')
max_egt_acc = get_accessions_from_prots(max_egt_prots, a_thaliana, 'Arabidopsis')
records, rec2locus, rec2organelle = get_records(max_egt_acc)

c_reinhardtii = parse_protein_map('../protein_maps_orthofinder/Chlamydomonas_reinhardtii.map')
qui_egt_prots = get_overlap(qui_gi, c_reinhardtii)
qui_egt_prots = format_prot_ids(qui_egt_prots, 'Chloroplastida')

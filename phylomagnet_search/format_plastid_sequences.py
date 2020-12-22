from Bio import SeqIO
import ete3
import glob
import re
import os

ncbi = ete3.ncbi_taxonomy.NCBITaxa()
tax2id = {'Prok': 2,
          'Colp29': 2584957,
          'Colp-29': 2584957,
          'Colp-38': 2583502,
          'Dinoflagellate': 2864,
          'Haptophyte': 2830,
          'Chrompodellid': 877183,
          'Chromosphaera': 127916}

for f in glob.glob("plastid_target_alns/*"):
    with open("plastid_target_faa/{}".format(os.path.basename(f)), 'w') as out:
        for rec in SeqIO.parse(f, 'fasta'):
            taxid = None
            for elem in reversed(re.split("-|_", rec.id)):
                try:
                    taxid = ncbi.get_name_translator([elem])[elem][0]
                except:
                    pass
                if taxid:
                    break
            if not taxid:
                for tax, tid in tax2id.items():
                    if tax in rec.id:
                        taxid = tid
            if not taxid:
                print(rec.id)
            rec.id = rec.id.replace('|', '_').replace('.', '_').replace('@', '_')
            rec.id = "{}.{}".format(taxid, rec.id)
            rec.seq = rec.seq.ungap('-')
            rec.description = ""
            SeqIO.write(rec, out, 'fasta')

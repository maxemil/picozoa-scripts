import ete3
import glob
from collections import defaultdict

ncbi = ete3.ncbi_taxonomy.NCBITaxa()

with open('contigs_NR.csv', 'w') as out:
    for f in glob.glob("*NR.out"):
        print(f)
        mito = defaultdict(bool)
        chloro = defaultdict(bool)
        taxids = defaultdict(set)
        for line in open(f):
            line = line.split('\t')
            if line[-1].strip():
                taxids[line[0]].update(line[-1].strip().split(';'))
            if 'mitochondri' in line[-2]:
                mito[line[0]] = True
            if 'plastid' in line[-2] or 'chloroplast' in line[-2]:
                chloro[line[0]] = True
        for contig, tids in taxids.items():
            failed_tids = []
            for tid in tids:
                try:
                    ncbi.get_taxid_translator([tid])[int(tid)]
                except:
                    failed_tids.append(tid)
            [tids.remove(t) for t in failed_tids]
            t = ncbi.get_topology(tids, intermediate_nodes=False)
            ancestor = ncbi.get_taxid_translator([t.name])[int(t.name)]
            mito_label = 'mitochondrial' if mito[contig] else '-'
            chloro_label = 'chloroplastic' if chloro[contig] else '-'
            print(f, contig, ancestor, mito_label, chloro_label, sep='\t', file=out)
            # print(f, contig, ancestor, mito_label, sep='\t')

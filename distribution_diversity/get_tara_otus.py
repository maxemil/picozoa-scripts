import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import seaborn as sns
# with open('Tara_OTUs.fasta', 'w') as out:
#     for line in open('TaraOceansV9_globaldataset2009-2012_otus_v20161202.tsv'):
#         if line.startswith('cid'):
#             continue
#         line = line.split('\t')
#         rec = SeqRecord(Seq(line[2]), id="OTU_{}".format(line[0]), description="")
#         SeqIO.write(rec, out, 'fasta')
# 
# 
# cmd = """/local/two/Software/vsearch-2.15.1/bin/vsearch --usearch_global \                    
         # Tara_OTUs.fasta --db Picozoa_18S_all_V9.fasta -iddef 1 --id 0.90 \
         # --blast6out TARA_Picozoa_90.out
# """
# subprocess.run(cmd.split())

picozoa_otus = [line.split()[0] for line in open("TARA_Picozoa_90.out")]

df = pd.read_csv('TaraOceansV9_globaldataset2009-2012_otus_v20161202.tsv', sep='\t')
df['OTU'] = df['cid'].apply(lambda x: "OTU_{}".format(x))
df = df[df['OTU'].apply(lambda x: x in picozoa_otus)]
df.index = df['OTU']
df = df.drop(['pid', 'refs', 'lineage', 'taxogroup','rtotab', 'totab', 'cid', 
             'md5sum', 'sequence', 'OTU'], axis=1)
pico_abundance = df.sum()

df_all = pd.read_csv('TaraOceansV9_globaldataset2009-2012_otus_v20161202.tsv', sep='\t')
df_all = df_all[df.columns]
all_abundance = df_all.sum()

env = pd.read_csv('datasets/TARA_ENV_DEPTH_SENSORS.tab', sep='\t', skiprows=2597)
env.index = env['Sample ID (TARA_barcode#, registered at ...)']
env = env[['Latitude', 'Longitude']]

def get_lat_lon(id, ll):
    try:
        return env.loc[id.split(',')[1], ll]
    except:
        print(id)
        return None

pico_rel = pico_abundance/all_abundance * 100
pico_rel = pd.DataFrame(pico_rel, columns=['Abundance'])
pico_rel['Latitude'] = pico_rel.index.map(lambda x: get_lat_lon(x, 'Latitude'))
pico_rel['Longitude'] = pico_rel.index.map(lambda x: get_lat_lon(x, 'Longitude'))
pico_rel = pico_rel.sort_values(by='Abundance')
pico_rel = pico_rel.drop_duplicates(subset=['Latitude','Longitude'],keep='last')
pico_rel = pico_rel[pico_rel['Abundance'] > 0.0]

def color_levels(abundance):
    if abundance < 1:
        return 'slategrey'
    elif abundance < 5:
        return 'darkviolet'
    elif abundance < 10:
        return 'midnightblue'
    elif abundance < 20:
        return 'orangered'
    elif abundance >= 20:
        return 'darkred'
pico_rel['color'] =  pico_rel['Abundance'].apply(lambda x: color_levels(x))


fig, ax = plt.subplots(figsize=(20,10), frameon=False)
m = Basemap(projection='cyl', llcrnrlat=-75, urcrnrlat=75,
            llcrnrlon=-165, urcrnrlon=85, resolution='l', ax=ax)
m.shadedrelief()
# m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
# m.drawmapboundary()
# m.drawcoastlines(color='gray', linestyle='dotted')
scat = m.scatter(pico_rel['Longitude'], pico_rel['Latitude'], latlon=True,
                c=pico_rel['Abundance'], cmap=sns.color_palette("light:r", as_cmap=True),
                alpha=0.5, s=pico_rel['Abundance']*15, norm=LogNorm())
                 # c=pico_rel['color'], alpha=0.5, s=pico_rel['Abundance']*15)
m.scatter(-122.013, 36.748, latlon=True, c='black', alpha=1, s=40, marker="^")
m.scatter(-122.357, 36.695, latlon=True, c='black', alpha=1, s=40, marker="^")
m.scatter(-123.49, 36.126, latlon=True, c='black', alpha=1, s=40, marker="^")
m.scatter(19.8633, 58.4880, latlon=True, c='black', alpha=1, s=40, marker="^")
# 'gist_heat_r' pico_rel['col'] sns.color_palette("Spectral_r", as_cmap=True)
# plt.colorbar(scat, label='relative abundance [%]', ax=ax)
for a, l in zip([0.1, 1, 10, 30], ["0.1%","1%","10%","30%"]):
    plt.scatter([], [], c='k', alpha=0.5, s=a*15,
                label=l, norm=LogNorm(), cmap=sns.color_palette("light:r", as_cmap=True))
plt.scatter([], [], c='black', alpha=1, s=40, label='SAG sampling', marker="^")
plt.legend(scatterpoints=1, frameon=False,
           labelspacing=1, loc='lower left')
plt.tight_layout()
plt.savefig('Tara_abundances.pdf')
plt.close()

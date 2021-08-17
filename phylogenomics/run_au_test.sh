#remove support values, e.g. in ete3:
# import ete3
# tree = ete3.PhyloTree("Picozoa_67taxa_317genes_BMGE_LGPMSF_BS.new", format=0) 
# tree.write(outfile='Picozoa_67taxa_317genes_BMGE_LGPMSF_BS.topology.tree', format=9)
# manually change topology
cat Picozoa_67taxa_317genes_BMGE_LGPMSF_BS.topology* > topologies2test.tree
sed  -i 's/\./_/g' topologies2test.tree
sed  -i "s/'//g" topologies2test.tree

cp 67_selection/67taxa_317genes_BMGE/concatenation/Picozoa_67taxa_317genes_BMGE.aln .
sed  -i 's/\./_/g' Picozoa_67taxa_317genes_BMGE.aln 
sed  -i 's/\//_/g' Picozoa_67taxa_317genes_BMGE.aln 

iqtree2 -s 67_selection/67taxa_317genes_BMGE/concatenation/Picozoa_67taxa_317genes_BMGE.aln \
       -z topologies2test.tree \
       -n 0 \
       -m LG+C60+F \
       -zb 10000 \
       -au \
       -zw \
       -nt 20 \
       -mem 70GB 
       # 
       # Tree      logL    deltaL  bp-RELL    p-KH     p-SH    p-WKH    p-WSH       c-ELW       p-AU
       # -------------------------------------------------------------------------------------------
       #   (Picozoa,(Rhodelphis,Rhodophyta)) -3576101.977       0   0.786 +  0.797 +      1 +  0.797 +  0.892 +     0.784 +    0.817 + 
       #   (Rhodelphis,(Picozoa,Rhodophyta))  -3576116.89  14.912   0.202 +  0.203 +   0.29 +  0.203 +  0.323 +     0.204 +    0.236 + 
       #   (Rhodophyta,(Rhodelphis,Picozoa)) -3576132.624  30.647  0.0114 - 0.0262 - 0.0467 - 0.0262 - 0.0569 +    0.0116 -   0.0362 - 
       # 
       # 
       # 

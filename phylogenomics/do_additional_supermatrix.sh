# only genes with 2 or more picozoa sequences in 67 dataset
for aln in 67_at_least_2_picozoa/*
do 
    if [ $(grep -c 'Picozoa' $aln) -lt 2 ]
    then 
        rm $aln
    fi
done

python ~software/misc-scripts/concatenate.py \
    -out 67_at_least_2_picozoa.aln \                         
    -s \
    -t 5 \
    -sep '@' \                  
    67_at_least_2_picozoa/*

iqtree2 -s 67_at_least_2_picozoa.aln -pre 67_at_least_2_picozoa \
        -bb 1000 -bnni -nt 20 -m LG+C60+F -mem 70G

python /local/two/Software/misc-scripts/color_tree.py \
        -i 67_at_least_2_picozoa.treefile \
        -o 67_at_least_2_picozoa.nex \
        -c ../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt


# fast site removal 
singularity exec ~software/phylofisher/phylofisher.sif fast_site_remover.py \
            -m Picozoa_67taxa_317genes_BMGE_fixnames.aln \
            -tr Picozoa_67taxa_317genes_BMGE_LGPMSF_BS.treefile \
            -of fasta -o fast_site_remover_5k -s 5000

#for each step until 50000 sites removed
for aln in step0.fas step1.fas step2.fas step3.fas step4.fas step5.fas step6.fas step7.fas step8.fas step9.fas step10.fas; 
do 
  iqtree2 -s $aln -pre ${aln%%.fas} -bb 1000 -nt 20 -m LG+C60+F;
done

for tree in *.treefile; 
do 
  python /local/two/Software/misc-scripts/color_tree.py \
          -i $tree \
          -o ${tree%%.treefile}.nex \
          -c ../../../../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt
done

# chi2-pruning

singularity exec ~software/alnpruner/alignment_pruner.sif alignment_pruner.pl \
              --file Picozoa_67taxa_317genes_BMGE.aln \
              --chi2_prune f0.25 > Picozoa_67taxa_317genes_BMGE_25chi2.aln

iqtree2 -s Picozoa_67taxa_317genes_BMGE_25chi2.aln -pre Picozoa_67taxa_317genes_BMGE_25chi2 \
        -bb 1000 -nt 20 -m LG+C60+F -mem 60G

python /local/two/Software/misc-scripts/color_tree.py \
        -i Picozoa_67taxa_317genes_BMGE_25chi2.treefile \
        -o Picozoa_67taxa_317genes_BMGE_25chi2.nex \
        -c ../../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt


singularity exec ~software/alnpruner/alignment_pruner.sif alignment_pruner.pl \
              --file Picozoa_67taxa_317genes_BMGE.aln \
              --chi2_prune f0.5 > Picozoa_67taxa_317genes_BMGE_50chi2.aln

iqtree2 -s Picozoa_67taxa_317genes_BMGE_50chi2.aln -pre Picozoa_67taxa_317genes_BMGE_50chi2 \
        -bb 1000 -nt 20 -m LG+C60+F -mem 50G
        
python /local/two/Software/misc-scripts/color_tree.py \
        -i Picozoa_67taxa_317genes_BMGE_50chi2.treefile \
        -o Picozoa_67taxa_317genes_BMGE_50chi2.nex \
        -c ../../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt
        

# Astral

for aln in *.mafft;
do
  iqtree2 -s $aln -pre ${aln%%.fas} \
        -bb 1000 -nt 4 -m TEST -mset LG -mrate G,R4 \
        -madd LG4X,LG4X+F,LG4M,LG4M+F -wbtl
done

for t in *.treefile; 
do
  python /local/two/Software/misc-scripts/color_tree.py \
          -i $t -o ${t%%.treefile}.nex \
          -c ../../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt
done


cat ../sgt_bootstraps/*.contree > sgt_all.new
sed -i 's/\//_/g' sgt_all.new
mkdir bootstraps
cp ../sgt_bootstraps/*.ufboot bootstraps
sed -i 's/\//_/g' bootstraps/*
for f in bootstraps/*.ufboot; do echo $f >> sgt_ufboot.txt; done

#clean rec names for astral:
"""
import glob
from Bio import SeqIO

recdict = {}
for f in glob.glob('../sgt_bootstraps/*.mafft'):
    for rec in SeqIO.parse(f, 'fasta'):
        recid = rec.id.replace('/', '_')
        recdict[recid.replace('@', '_')] = recid.split('@')[0]

with open("replace_names.txt", 'w') as out:
    for k, v in recdict.items():
        print(k, v, sep='\t', file=out)
"""
python /local/two/Software/misc-scripts/replace_names.py  \
          -i replace_names.txt \
          -f sgt_all.new -b

parallel -j 15 python /local/two/Software/misc-scripts/replace_names.py \
          -i replace_names.txt -f {} ::: bootstraps/*.ufboot

singularity exec ~software/phylofisher/phylofisher.sif \
          java -Xmx40000M -jar /opt/conda/envs/phylofisher/ASTRAL-5.7.3/astral.5.7.3.jar \
          -i sgt_all.new \
          -b sgt_ufboot.txt \
          -o astral_tree.new 2> astral_tree.log

tail -n 1 astral_tree.new >astral_tree.main_tree

python /local/two/Software/misc-scripts/color_tree.py \
        -i astral_tree.main_tree \
        -o astral_tree.nex \
        -c ../../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt

singularity exec ~software/phylofisher/phylofisher.sif \
          java -Xmx40000M -jar /opt/conda/envs/phylofisher/ASTRAL-5.7.3/astral.5.7.3.jar \
          -i sgt_all.new \
          -q astral_tree.main_tree \
          -o astral_tree_pp.new 2>astral_tree_pp.log
        
        
# run modeltest on all new alignments
for f in *.aln; 
do
  /opt/iqtree-2.1.2-Linux/bin/iqtree2 \
        -s $f \
        -pre "${f%%.aln}"_Modeltest \
        -nt 20 \
        -m TESTONLY \
        -mset LG \
        -madd LG4X+F,LG+C20+F,LG+C40+F,LG+C60+F,LG4X,LG+C20,LG+C40,LG+C60 \
        -mrate I,G,G+I \
        -mem 100G
done


# remove genes with just a single picozoa
mkdir 67_multiple_picozoa
cp ../25_phylogenomics_taxon_sampling/67_selection/67taxa_317genes_BMGE/alignments/* 67_multiple_picozoa
rm 67_multiple_picozoa/*.mafft

for aln in 67_multiple_picozoa/*
do 
    if [ $(grep -c 'Picozoa' $aln) -eq 1 ]
    then 
        rm $aln
    fi
done

cd 67_multiple_picozoa

python ~software/misc-scripts/concatenate.py \
    -out 67_multiple_picozoa.aln \
    -s \
    -t 5 \
    -sep '@' *.aln

iqtree2 -s 67_multiple_picozoa.aln -pre 67_multiple_picozoa \
        -bb 1000 -nt 20 -m LG+C60+F -mem 70G

python /local/two/Software/misc-scripts/color_tree.py \
        -i 67_multiple_picozoa.treefile \
        -o 67_multiple_picozoa.nex \
        -c ../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt


# remove genes with no picozoa
mkdir 67_only_picozoa
cp ../25_phylogenomics_taxon_sampling/67_selection/67taxa_317genes_BMGE/alignments/* 67_only_picozoa
rm 67_only_picozoa/*.mafft

for aln in 67_only_picozoa/*
do 
    if [ $(grep -c 'Picozoa' $aln) -eq 0 ]
    then 
        rm $aln
    fi
done

cd 67_only_picozoa

python ~software/misc-scripts/concatenate.py \
    -out 67_only_picozoa.aln \
    -s \
    -t 5 \
    -sep '@' *.aln

iqtree2 -s 67_only_picozoa.aln -pre 67_only_picozoa \
        -bb 1000 -nt 20 -m LG+C60+F -mem 40G

python /local/two/Software/misc-scripts/color_tree.py \
        -i 67_only_picozoa.treefile \
        -o 67_only_picozoa.nex \
        -c ../25_phylogenomics_taxon_sampling/67_selection/taxa_colors.txt

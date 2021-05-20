nextflow run align_tree_clusters.nf \
                    --faa "Orthogroup_Sequences/OG0000[3-6]*.fa" \
                    -with-singularity gene_trees.sif \
                    -resume
for f in *.treefile; 
do 
  python rename_orthogroup_tree.py $f Picozoa taxonomy; 
done

for f in *.Picozoa.treefile; 
do 
  python /local/two/Software/misc-scripts/color_tree.py -i $f \
                            -o ${f%%.treefile}.nex -c taxa_colors_Picozoa.txt; 
done

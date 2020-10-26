taxon="Glaucophyta"

mkdir Orthogroup_Selection_$taxon/trees
mkdir Orthogroup_Selection_$taxon/orthogroups

mv Orthogroup_Selection_$taxon/OG*.fa Orthogroup_Selection_$taxon/orthogroups

for fa in Orthogroup_Selection_$taxon/orthogroups/*;
do
  og=$(basename ${fa%%.fa})
  cp Orthogroup_alignments_trees/trees/$og.treefile Orthogroup_Selection_$taxon/trees/
done

cp taxa_colors_orthogroups.txt Orthogroup_Selection_$taxon/taxa_colors_$taxon.txt
sed -i '/^Picozoa/ s/red/darkred/g' Orthogroup_Selection_$taxon/taxa_colors_$taxon.txt
echo "$taxon\t red" >> Orthogroup_Selection_$taxon/taxa_colors_$taxon.txt

parallel -j 40 python rename_orthogroup_tree.py {} Helicosporidium ::: Orthogroup_Selection_$taxon/trees/*[0-9].treefile
parallel --link -j 40 python /local/two/Software/misc-scripts/color_tree.py -i {2} -o {1} -c Orthogroup_Selection_$taxon/taxa_colors_$taxon.txt ::: $(ls Orthogroup_Selection_$taxon/trees/*.$taxon.treefile | sed 's/.treefile/.nex/g') ::: $(ls Orthogroup_Selection_$taxon/trees/*.$taxon.treefile)

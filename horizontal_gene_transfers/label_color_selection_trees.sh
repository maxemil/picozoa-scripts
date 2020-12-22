for dir in Orthogroup_LGT_selections/* ;
do
  taxon=$(echo $dir | sed 's/Orthogroup_LGT_selections\///')
  echo $taxon
  # mkdir $dir/fast_trees
  # mkdir $dir/orthogroups
  # 
  # mv $dir/OG*.fa $dir/orthogroups
  
  for fa in $dir/orthogroups/*;
  do
    og=$(basename ${fa%%.fa})
    cp Orthogroup_LGT_fast_trees/trees/$og.treefile $dir/fast_trees/
  done

  cp taxa_colors_orthogroups.txt $dir/taxa_colors_$taxon.txt
  sed -i '/^Picozoa/ s/red/darkred/g' $dir/taxa_colors_$taxon.txt
  echo "$taxon\t red" >> $dir/taxa_colors_$taxon.txt
  if [ $taxon = 'Picozoa' ];
  then
    parallel -j 40 python rename_orthogroup_tree.py {} $taxon taxonomy ::: $dir/fast_trees/*[0-9].treefile
  else
    parallel -j 40 python rename_orthogroup_tree.py {} $taxon Name ::: $dir/fast_trees/*[0-9].treefile
  fi
  parallel --link -j 40 python /local/two/Software/misc-scripts/color_tree.py -i {2} -o {1} -c $dir/taxa_colors_$taxon.txt ::: $(ls $dir/fast_trees/*.$taxon.treefile | sed 's/.treefile/.nex/g') ::: $(ls $dir/fast_trees/*.$taxon.treefile)
done

sed -i 's/Picozoa..SAG33/Cryptophyceae..SAG33/g' Orthogroup_LGT_selections/Picozoa/fast_trees/*.Picozoa.*

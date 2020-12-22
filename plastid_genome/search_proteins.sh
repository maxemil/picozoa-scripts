wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.protein.faa.gz
diamond makedb --in plastid.3.protein.faa --db plastid.3.protein.dmnd --threads 40

for f in ../../21_contamination_final_assemblies/*.clean.faa; 
do
  base=$(basename $f .clean.faa); 
  diamond blastp --query $f --db plastid.3.protein.dmnd --out $base.out \
                  --threads 20 --outfmt 6 --more-sensitive --evalue 1e-5; 
done



for assembly in ../../09_coassemblies/*/spades/*.gfa; 
do
  cosag=${assembly%%/spades*} 
  cosag=${cosag##*/}
  singularity exec -B /local getorganelle.sif get_organelle_from_assembly.py \
                                -F embplant_pt,other_pt \
                                -g $assembly \
                                -o $cosag"_getOrganelle_ass" \
                                -t 4
done

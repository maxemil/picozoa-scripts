proteins="Rhodelphis_Rhodophytes.faa"

mkdir added_target_faa
mkdir added_target_hmm
for f in references_eggnog/*/*.aln
do 
  base=$(basename ${f%%.*})
  hmm="added_target_hmm/$base.hmm"
  hmmbuild -n $base $hmm $f > /dev/null
  hmmsearch  --tblout added_target_hmm/$base.out -E 1e-05 $hmm $proteins > /dev/null
  python extract_hmm_hits.py added_target_hmm/$base.out $proteins added_target_faa/$base.fasta
  cat eggnog_target_faa/$base.fasta >> added_target_faa/$base.fasta  
done

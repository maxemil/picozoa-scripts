cat pr2_picozoa_99.fasta pr2_katablepharids.90.fasta pr2_Cryptophyta.fasta Yoon_18S.fasta 18S_picozoa_SAGs_euk.fasta > PK.fasta
sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picozoa_XXXX|Picozoa_XXXX_sp.//g' PK.fasta
sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picomonas|//g' PK.fasta
sed -i 's/|18S_rRNA|nucleus|/_/g' PK.fasta

mafft-einsi --thread 20 --adjustdirection PK.fasta > PK.aln
sed -i 's/_R_//g' PK.aln

iqtree -s PK.aln -pre PK -bb 1000 -bnni -nt 10 -m GTR+G -redo

# 
# cat pr2_picozoa_99.fasta pr2_katablepharids.90.fasta pr2_Cryptophyta.fasta Yoon_18S.fasta initial_SAGs_18S.fasta > initial_SAGs_ref.fasta
# sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picozoa_XXXX|Picozoa_XXXX_sp.//g' initial_SAGs_ref.fasta
# sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picomonas|//g' initial_SAGs_ref.fasta
# sed -i 's/|18S_rRNA|nucleus|/_/g' initial_SAGs_ref.fasta
# 
# # bfg -v --search-records "HQ865994.1.881_U_clone_SI021806_215_SGSI941" -f initial_SAGs_ref.fasta > initial_SAGs_ref_clean.fasta
# 
# mafft-einsi --thread 20 --adjustdirection initial_SAGs_ref.fasta > initial_SAGs_ref.aln
# sed -i 's/_R_//g' initial_SAGs_ref.aln
# 
# iqtree -s initial_SAGs_ref.aln -pre initial_SAGs_ref -bb 1000 -bnni -nt 10 -m GTR+G -redo

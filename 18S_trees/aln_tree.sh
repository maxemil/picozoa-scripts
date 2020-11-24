cat pr2_picozoa_99.fasta pr2_katablepharids.90.fasta pr2_Cryptophyta.fasta Yoon_18S.fasta 18S_picozoa_SAGs_euk.fasta > PK.fasta
sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picozoa_XXXX|Picozoa_XXXX_sp.//g' PK.fasta
sed -i 's/|Eukaryota|Hacrobia|Picozoa|Picozoa_X|Picozoa_XX|Picozoa_XXX|Picomonas|//g' PK.fasta
sed -i 's/|18S_rRNA|nucleus|/_/g' PK.fasta

mafft-einsi --thread 20 --adjustdirection PK.fasta > PK.aln
sed -i 's/_R_//g' PK.aln

iqtree -s PK.aln -pre PK -bb 1000 -bnni -nt 10 -m GTR+G -redo

/opt/diamond-2.0.6/diamond makedb --in 210113_nr \
                                  --db 210113_nr.dmnd \
                                  --threads 4 \
                                  --taxonmap prot.accession2taxid

for f in *.faa
do
   base=${f%%.faa}
   echo $base
   /opt/diamond-2.0.6/diamond blastp -d /media/Data_2/nr_fasta/210113_nr.dmnd \
           -q $f \
           -o $base"_NR.out" \
           --more-sensitive \
           --quiet \
           -e 0.001 \
           --top 10 \
           --threads 4 \
           --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids
done

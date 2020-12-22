#!/bin/bash -l

#SBATCH -A snic2020-15-58 
#SBATCH -p core 
#SBATCH -n 15
#SBATCH -t 48:00:00
#SBATCH -J diamond-nr

module load bioinfo-tools
module load diamond/2.0.4

for faa in ../12_final_assemblies/*.faa; 
do
        diamond blastp -q $faa \
                      --db $DIAMOND_NR \
                      -e 1e-5 \
                      --out "$faa".out \
                      --outfmt 6 qseqid sseqid staxids evalue pident length mismatch gapopen qstart qend sstart send bitscore \
                      --threads 15 \
                      --top 2
done

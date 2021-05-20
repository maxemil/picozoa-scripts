BURNIN=1500
CYC=3600

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c1 \
        *chain1.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c2 \
        *chain2.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c3 \
        *chain3.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c1_c2 \
        *chain1.treelist *chain2.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c1_c3 \
        *chain1.treelist *chain3.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_c2_c3 \
        *chain2.treelist *chain3.treelist

singularity exec -B $PWD PhyloBayesMPI.sif bpcomp -x $BURNIN 1 $CYC -o "$CYC"cyc_all \
        *chain1.treelist *chain2.treelist *chain3.treelist

for t in $CYC*.con.tre; do 
        tree=${t%.tre}".nex"
        python /local/two/Software/misc-scripts/color_tree.py \
                    -i $t \
                    -o $tree \
                    -c ../../taxa_colors.txt
done

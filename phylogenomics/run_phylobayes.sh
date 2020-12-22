ALIGNMENT=$1
BASE=$(basename ${ALIGNMENT%%.*})
trimal -in $ALIGNMENT  -out $BASE.phy -phylip

echo $BASE
singularity exec -B $PWD PhyloBayesMPI.sif mpirun -np 8 pb_mpi \
                                      -d $BASE.phy \
                                      -cat -gtr \
                                      -dc $BASE"_PB_chain1" &

singularity exec -B $PWD PhyloBayesMPI.sif mpirun -np 8 pb_mpi \
                                      -d $BASE.phy \
                                      -cat -gtr \
                                      -dc $BASE"_PB_chain2" &

singularity exec -B $PWD PhyloBayesMPI.sif mpirun -np 8 pb_mpi \
                                      -d $BASE.phy \
                                      -cat -gtr \
                                      -dc $BASE"_PB_chain3" &

#! /bin/bash

#declare -a nprocs=(1 2 4 8 16 32 64)

#echo ${nprocs[@]}

#for nprocs in 1 2 4 8 16 32 64
for nprocs in 8
do
    echo -e "64\n $nprocs \n" | cat >proc.input
    ./proc <proc.input
    for preconditioner in asm hypre ml
    do
        numactl --localalloc --cpunodebind=4,5,6,7 mpirun -n $nprocs ./solver -pc_type $preconditioner -log_summary logfiles_numericalEfficiency_$preconditioner/log.64.$nprocs
        cp ERR.out logfiles_numericalEfficiency_$preconditioner/ERR.64.$nprocs
    done
done

echo DONE!

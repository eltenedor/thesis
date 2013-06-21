#! /bin/bash

mkdir -p logfiles_parallelEfficiency_mpich_bjacobi
mkdir -p logfiles_parallelEfficiency_mpich_asm
mkdir -p logfiles_parallelEfficiency_mpich_hypre
mkdir -p logfiles_parallelEfficiency_mpich_ml

preconditioner[0]=bjacobi
preconditioner[1]=asm
preconditioner[2]=hypre
preconditioner[3]=ml

its[0]=193
its[1]=192
its[2]=6
its[3]=8

for nprocs in 1 2 4 8 16 32 64
do
    	echo -e "64\n $nprocs \n" | cat >proc.input
    	./proc <proc.input
	for rankfilename in ~/userfiles_mpich/rankfile.${nprocs}.*
	do
		user=$(<${rankfilename})
		for i in 0 1 2 3 
		do 
			echo -e "mpiexec.hydra -np ${nprocs} --binding user:${user} ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_parallelEfficiency_mpich_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}"
			mpiexec.hydra -np ${nprocs} --binding user:${user} ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_parallelEfficiency_mpich_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}
			cp ERR.out logfiles_parallelEfficiency_mpich_${preconditioner[${i}]}/ERR.128.${rankfilename#*.}
		done
	done
done

echo DONE!

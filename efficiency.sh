#! /bin/bash

mkdir -p logfiles_parallelEfficiency_bjacobi
mkdir -p logfiles_parallelEfficiency_asm
mkdir -p logfiles_parallelEfficiency_hypre
mkdir -p logfiles_parallelEfficiency_ml

preconditioner[0]=bjacobi
preconditioner[1]=asm
preconditioner[2]=hypre
preconditioner[3]=ml

its[0]=672
its[1]=672
its[2]=6
its[3]=11

for nprocs in 1 2 4 8 16 32 64
do
    	echo -e "64\n $nprocs \n" | cat >proc.input
    	./proc <proc.input
	for rankfilename in ~/rankfiles/rankfile.${nprocs}.*
	do
		echo $rankfilename
		cat ${rankfilename}
		for i in 0 1 2 3 
		do 
			echo -e "mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_parallelEfficiency_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}"
			mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_parallelEfficiency_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}
			cp ERR.out logfiles_parallelEfficiency_${preconditioner[${i}]}/ERR.128.${rankfilename#*.}
		done
	done
done

echo DONE!

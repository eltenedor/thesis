#! /bin/bash

mkdir -p logfiles_scaling_bjacobi
mkdir -p logfiles_scaling_ml

preconditioner[0]=bjacobi
preconditioner[1]=ml

	for rankfilename in ~/rankfiles/rankfile.1.a # ~/rankfiles/rankfile.2.b ~/rankfiles/rankfile.4.c ~/rankfiles/rankfile.8.c ~/rankfiles/rankfile.16.d ~/rankfiles/rankfile.32.c
	do
		nprocs=${rankfilename#*.}
		nprocs=${nprocs%.*}
		echo -e "64\n $nprocs \n" | cat >proc.input
		./proc <proc.input
		#echo nprocs
		#echo $nprocs
		#echo $rankfilename
		#cat ${rankfilename}

		for i in 0 # 1 
		do 
			echo -e "mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_scaling_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}"
			mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_scaling_${preconditioner[${i}]}/log.128.${rankfilename#*.} -ksp_max_it ${its[${i}]}
			cp ERR.out logfiles_scaling_${preconditioner[${i}]}/ERR.128.${rankfilename#*.}
		done
	done

echo DONE!

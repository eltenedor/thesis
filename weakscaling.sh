#! /bin/bash

mkdir -p logfiles_weakscaling_bjacobi
mkdir -p logfiles_weakscaling_ml

preconditioner[0]=bjacobi
preconditioner[1]=ml

	for rankfilename in ~/rankfiles/rankfile.1.a ~/rankfiles/rankfile.2.b ~/rankfiles/rankfile.4.c ~/rankfiles/rankfile.8.c ~/rankfile/rankfile.16.d
	do
		nprocs=${rankfilename#*.}
		nprocs=${nprocs%.*}
		echo -e "$nprocs\n $nprocs \n" | cat >proc.input
		./proc <proc.input
		cp weakScalingInput/gridfiles/${nprocs}B/* .
		echo -e "$nprocs\n grid \n" | cat >gr.input
		numactl --localalloc --cpunodebind=4 ./grgen <gr.input
		cd ../pet_src
		make clean all -f premakefile
		cd ../128_weakScaling
		echo -e "$nprocs\n" | cat >pre.input
		echo Preprocessing
		numactl --localalloc --cpunodebind=4 ./preprocessing <pre.input
		cd ../pet_src
		make clean all -f somakefile
		cd ../128_weakScaling

		#echo nprocs
		#echo $nprocs
		#echo $rankfilename
		#cat ${rankfilename}

		for i in 0 1 
		do 
			echo -e "mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_weakscaling_${preconditioner[${i}]}/log.128.${rankfilename#*.} "
			mpirun -np ${nprocs} -rf $rankfilename ./solver -pc_type ${preconditioner[${i}]} -log_summary logfiles_weakscaling_${preconditioner[${i}]}/log.128.${rankfilename#*.} 
			cp ERR.out logfiles_weakscaling_${preconditioner[${i}]}/ERR.128.${rankfilename#*.}
		done
	done

echo DONE!

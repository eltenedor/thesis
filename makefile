all: solver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# default values, will be overwritten by bash
expath = ../order_verification/blockStructured/test/

solver_src = \
	${expath}parameterModule.F90 \
	controlModule.F90 \
	boundaryModule.F90 \
	characterModule.F90 \
	fluxModule.F90 \
	geoModule.F90 \
	gradModule.F90 \
	indexModule.F90 \
	varModule.F90 \
	scModule.F90 \
	mmsModule.F90 \
	petsc_ksp_module.F90 \
	coefModule.F90 \
	solver.F90 \
    
	
solver_obj = ${solver_src:%.F90=%.o}
solver_mod = ${solver_src:%.F90=%.mod}

# Executables

solver: ${solver_obj} chkopts
	-${FLINKER} -o ${expath}/solver ${solver_obj} ${PETSC_LIB}
	${RM} -f ${solver_obj} ${solver_mod}
	
cln:
	${RM} -f *.o *.mod $(expath)*.o $(expath)*.mod
	
	

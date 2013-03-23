all: solver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

# default values, will be overwritten by bash
resolution = 2
expath = ../order_verification/blockStructured/${resolution}

solver_src = \
	${expath}/paramMod.F90 \
	logicMod.F90 \
	bcMod.F90 \
	chMod.F90 \
	fluxMod.F90 \
	geoMod.F90 \
	gradMod.F90 \
	indMod.F90 \
	preProcIndMod.F90 \
	varMod.F90 \
	scMod.F90 \
	mmsMod.F90 \
	petsc_ksp_module.F90 \
	coefMod.F90 \
	solver.F90 \
    
grgen_src = \
	paramIng.F90 \
	bcMod.F90 \
	chMod.F90 \
	geoMod.F90 \
	indMod.F90 \
	grgen.F90 \

preprocessing_src = \
	${expath}/paramMod.F90 \
	bcMod.F90 \
	geoMod.F90 \
	preProcIndMod.F90 \
	preprocessing.F90 \
	
solver_obj = ${solver_src:%.F90=%.o}
solver_mod = ${solver_src:%.F90=%.mod}
grgen_obj = ${grgen_src:%.F90=%.o}
preprocessing_obj = ${preprocessing_src:%.F90=%.o}
preprocessing_mod = ${preprocessing_src:%.F90=%.mod}

# Executables

solver: ${solver_obj} chkopts
	-${FLINKER} -o ${expath}/solver ${solver_obj} ${PETSC_LIB}
	${RM} -f ${solver_obj} ${solver_mod}
	
grgen: ${grgen_obj} chkopts
	gfortran -o ${expath}/grgen ${grgen_obj}
	${RM} -f *.o *.mod
	
preprocessing: ${preprocessing_obj} chkopts
	gfortran -o ${expath}/preprocessing ${preprocessing_obj}
	${RM} -f ${preprocessing_obj} ${preprocessing_mod}
	
cln:
	${RM} -f *.o *.mod ${solver_obj} ${solver_mod} ${preprocessing_obj} ${preprocessing_mod}
	
	

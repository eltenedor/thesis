all: solver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

resolution = 2
expath = ../order_verification/petsc_cds/${resolution}

solver_src = \
	${expath}/paramMod.F90 \
    bcMod.F90 \
    chMod.F90 \
    fluxMod.F90 \
    geoMod.F90 \
    gradMod.F90 \
    indMod.F90 \
    varMod.F90 \
    scMod.F90 \
    mmsMod.F90 \
    logicMod.F90 \
    coefMod.F90 \
    solver.F90 \
    
grgen_src = \
	paramIng.F90 \
	bcMod.F90 \
	chMod.F90 \
	geoMod.F90 \
	indMod.F90 \
	grgen.F90 \
	
solver_obj = $(solver_src:%.F90=%.o)
grgen_obj = $(grgen_src:%.F90=%.o)


# Executables

solver: ${solver_obj} chkopts
	-${FLINKER} -o ${expath}/solver ${solver_obj} ${PETSC_LIB}
	${RM} ${solver_obj}
	
grgen: ${grgen_obj} chkopts
	gfortran -o ${expath}/grgen ${grgen_obj}
	${RM} ${grgen_obj}
	
cln:
	${RM} ${solver_obj} ${grgen_obj}
	
	

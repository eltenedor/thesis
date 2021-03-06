module controlModule

    implicit none
#include <finclude/petscsysdef.h>

    logical :: LTIME=.false. ! true for instationary problems
    logical :: CONVERGED=.false. ! stopping criterion in petsc_ksp_module
    logical :: STARTED,FOUND ! used for optimized preprocessing
    !integer :: pcInd, kspInd
    integer :: reasonInt, itsInt ! store convergedReason and numer of innerIterations
    integer,parameter :: PROCOFFSET=10,BLOCKOFFSET=20 ! used to determine inputfile names
    real*8 :: tges ! cpu time used for linear solve

    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    PetscScalar :: res_Scalar
    PetscLogDouble :: time1, time2

end module controlModule

#define PETSC_AVOID_DECLARATIONS
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#undef PETSC_AVOID_DECLARATIONS

!#########################################################
module petsc_ksp_module
!#########################################################

    implicit none

    KSP :: ksp
    PC :: pc
    PetscErrorCode :: ierr
    KSPConvergedReason :: reason
    Mat :: A2
    Vec :: vt1,vt2,res

contains

!#########################################################
subroutine setUpKSP
!#########################################################

    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

    ! Set runtime options

    call KSPSetFromOptions(ksp,ierr)

end subroutine setUpKSP

!#########################################################
subroutine solveSys(A,b,x,N,LS,tol)
!#########################################################

    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>

    Mat, intent(in) :: A
    Vec, intent(in) :: b
    Vec, intent(in out) :: x
    integer, intent(in) :: N,LS
    PetscScalar, intent(in out) :: tol


    if(LS.eq.1) then
        call MatConvert(A,MATSAME,MAT_INITIAL_MATRIX,A2,ierr)
        call VecDuplicate(x,vt1,ierr)
        call VecDuplicate(x,vt2,ierr)
        call VecDuplicate(x,res,ierr)
    else
        call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
    endif

    if(N.le.10) call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

    ! Set operators

    call KSPSetOperators(ksp,A,A2,SAME_PRECONDITIONER,ierr)

    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION, &
            & PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
    ! Solve the linear system

    call KSPSolve(ksp,b,x,ierr)
    call KSPGetConvergedReason(ksp,reason,ierr)
    !call PetscPrintf(PETSC_COMM_WORLD,reason,ierr);

    !call KSPBuildResidual(ksp,vt1,res,vres,ierr)
    call KSPInitialResidual(ksp,x,vt1,vt2,res,b,ierr)

    call VecMin(vt2,PETSC_NULL_INTEGER,tol,ierr)
    tol=abs(tol)

    !print *, TOL

    ! View solver info

    call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if(N.le.10) then
        call VecView(vt2,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    end if

end subroutine solveSys
    
!#########################################################
subroutine cleanUp(A,b,x)
!#########################################################

    
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>

    Mat, intent(in out) :: A
    Vec, intent(in out) :: b,x

    call VecDestroy(x,ierr)
    call VecDestroy(b,ierr)
    call VecDestroy(vt1,ierr)
    call VecDestroy(vt2,ierr)
    call VecDestroy(res,ierr)
    call MatDestroy(A,ierr)
    call MatDestroy(A2,ierr)
    call KSPDestroy(ksp,ierr) 

end subroutine cleanUp

end module petsc_ksp_module

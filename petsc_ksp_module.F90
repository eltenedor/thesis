!#########################################################
module petsc_ksp_module
!#########################################################

    implicit none
#include <finclude/petscsysdef.h>
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>
#include <finclude/petscpcdef.h>
#include <finclude/petsckspdef.h>

    KSP :: ksp
    PC :: pc
    PetscErrorCode :: ierr
    KSPConvergedReason :: reason
    PetscInt :: its
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

    use logic
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    Mat, intent(in) :: A
    Vec, intent(in) :: b
    Vec, intent(in out) :: x
    integer, intent(in) :: N,LS
    PetscScalar, intent(in out) :: tol
    PetscReal :: r2, b2, rtol, rinit
    
    ! Calculate initial Residual
    
    if(LS.eq.1) then
        ! Navier Stokes, muss die Matrix aktualisiert werden
        call MatConvert(A,MATSAME,MAT_INITIAL_MATRIX,A2,ierr)
        call VecDuplicate(x,vt1,ierr)
        call VecDuplicate(x,vt2,ierr)
        call VecDuplicate(x,res,ierr)
        !
        ! Set starting tolerance 1e-4
        !
        rtol = 1e-4
        !rtol = 1e-10
    else
        call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
        !
        ! Check convergence and calculate new relative Tolerance
        !
        call KSPInitialResidual(ksp,x,vt1,vt2,res,b,ierr)
        call VecMin(res,PETSC_NULL_INTEGER,tol,ierr)
        tol=abs(tol)
        if (tol<1e-12 .and. rtol < 1e-9) then
            CONVERGED=.true.
            print *, "Final tolerance: ", tol
            return
        else
            call VecNorm(vt2,NORM_2,r2,ierr)
            call VecNorm(b,NORM_2,b2,ierr)
            rtol=(r2/b2)/100.0            
        endif
    endif

    ! Set operators

    call KSPSetOperators(ksp,A,A2,SAME_PRECONDITIONER,ierr)

    call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION, &
            !& PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
            & PETSC_DEFAULT_DOUBLE_PRECISION,20000,ierr)

            
    ! Solve the linear system

    call KSPSolve(ksp,b,x,ierr)
    
    ! Get KSP information
    call KSPGetConvergedReason(ksp,reason,ierr)
    call KSPGetIterationNumber(ksp,its,ierr)

    reasonInt=reason
    itsInt=its

    !print *, TOL

    ! View solver info

    call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if(N.le.64) then
        print *, 'Matrix A:'
        call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !print *, 'Vector vt2:'
        !call VecView(vt2,PETSC_VIEWER_STDOUT_WORLD,ierr)
        !print *, 'Vector res:'
        !call VecView(res,PETSC_VIEWER_STDOUT_WORLD,ierr)
        print *, 'Vector x:'
        call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
        print *, 'Vector b:'
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
#include <finclude/petscpc.h>

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

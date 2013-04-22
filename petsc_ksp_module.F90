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
    KSPConvergedReason :: reason
    PetscInt :: its
    Mat :: A2
    Vec :: vt1,vt2,res
    PetscReal :: b_real, rtol, rinit, rfinal

    PetscErrorCode :: ierr

    Private::ierr

contains

!#########################################################
subroutine setUpKSP
!#########################################################

    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscksp.h>

    ! Create KSPContext

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

    ! Set runtime options

    call KSPSetFromOptions(ksp,ierr)


end subroutine setUpKSP

!=========================================================
!>  solves a linear system. Linear System Matrix A can be
!>  updated for each solve, but the same preconditioning 
!>  matrix is retained for all calls of this subroutine
!#########################################################
subroutine solveSys(A,b,x,N,LS,r_scalar)
!#########################################################

    use controlModule
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
    PetscScalar, intent(in out) :: r_scalar
    
    if(LS.eq.1) then
        ! save initial stiffness matrix inside the scope of the module
        ! to retain matrix between successive subroutine calls
        call MatConvert(A,MATSAME,MAT_INITIAL_MATRIX,A2,ierr)

        ! Initialize temporary work vectors
        call VecDuplicate(x,vt1,ierr)
        call VecDuplicate(x,vt2,ierr)
        call VecDuplicate(x,res,ierr)

        ! Calculate initial Residual
        call VecNorm(b,NORM_2,b_real,ierr)
        ! set final tolerance
        rfinal = b_real*1e-8
        rtol = 1e-2
        !rtol = 1e-8
        if (rank .eq. 0) print *, 'SETTING FINAL RESIDUAL TO: ', rfinal
        if (rank .eq. 0) print *, 'SETTING RELATIVE TOLERANCE TO: ', rtol
        if (rank .eq. 0) print *, 'INITIAL RESIDUAL: ', b_real
        r_scalar = b_real
    else
        call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
        !
        ! Check convergence by calculating the 2-Norm of the residual vector
        ! and calculate new relative Tolerance (
        ! (available are 1-,2- and INFTY-Norm; if the spectral radius of A < 1
        ! every Norm should converge to zero)
        !
        ! Calculate initial Residual
        call KSPInitialResidual(ksp,x,vt1,vt2,res,b,ierr)
        call VecNorm(vt2,NORM_2,rinit,ierr)

        if (rinit<rfinal) then
            CONVERGED=.true.
            if (rank .eq. 0) print *, 'FINAL RESIDUAL: ', rinit
            if (rank .eq. 0) write (9,*) 'FINAL RESIDUAL: ', rinit
            r_scalar = rinit
            return
        else
            call VecNorm(b,NORM_2,b_real,ierr)
            rtol=(rinit/b_real)*1e-2
            if (rank .eq. 0) print *, 'MOMENTARY RESIDUAL: ', rinit
            !if (rank .eq. 0) print *, "SETTING RELATIVE TOLERANCE TO: ", rtol
            r_scalar = rinit
        endif
    endif

    ! Set operators

    call KSPSetOperators(ksp,A,A2,SAME_PRECONDITIONER,ierr)

    call KSPSetTolerances(ksp,rtol,PETSC_DEFAULT_DOUBLE_PRECISION, &
            !& PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
            & PETSC_DEFAULT_DOUBLE_PRECISION,20000,ierr)
            
    ! Solve the linear system

    call PetscGetTime(time1,ierr)
    call KSPSolve(ksp,b,x,ierr)
    call PetscGetTime(time2,ierr)
    
    ! Get KSP information
    call KSPGetConvergedReason(ksp,reason,ierr)
    call KSPGetIterationNumber(ksp,its,ierr)

    reasonInt=reason
    itsInt=its

    if (rank .eq. 0) print '(A,I6,A,I2)', '... REACHED AFTER: ', its, ', REASON: #', reason

    !print *, TOL
    ! View solver info
    !call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if (N.le.16) then
        call PetscObjectSetName(A,'Matrix A:',ierr)
        call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call PetscObjectSetName(vt2,'Vector vt2:',ierr)
        call VecView(vt2,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call PetscObjectSetName(res,'Vector res:',ierr)
        call VecView(res,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call PetscObjectSetName(x,'Vector x:',ierr)
        call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call PetscObjectSetName(b,'Vector b:',ierr)
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

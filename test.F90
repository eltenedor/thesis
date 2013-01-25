program main

    implicit none
    
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

    real*8 :: R1,R2
    integer,parameter :: N=4
    Mat :: A
    PetscInt :: col(4),row(4),one,four
    PetscScalar :: val(4),Q,AP,AN,AW,AE
    PetscErrorCode :: ierr
    
    
    
    four=4
    one=1
    
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    
    call MatCreate(PETSC_NULL_CHARACTER,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
    call MatSetFromOptions(A,ierr)
    call MatSetUp(A,ierr)
    
    AE=1.0
    
        col(1)=0
        col(2)=1
        col(3)=1
        col(4)=2
        row(1)=0
        row(2)=0
        row(3)=1
        row(4)=1
        val(1)=1.0
        val(2)=1.0
        val(3)=1.0
        val(4)=1.0
        
        !call MatSetValues(A,four,row,four,col,val,INSERT_VALUES,ierr)
        
        !call MatSetValue(A,one,one,AE,INSERT_VALUES,ierr)
        
        !call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
        
        !call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
        
        !call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
        
    call PetscFinalize(ierr)
print *, 'START'
end program main

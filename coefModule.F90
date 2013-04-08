module coefModule

    use parameterModule
    implicit none
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>


    real(KIND=PREC) :: Q(NXYZA), AP(NXYZA), AS(NXYZA), AN(NXYZA),&
                       AW(NXYZA), AE(NXYZA), AB(NXYZA), AT(NXYZA),&
                       AF(NFACEAL)
                       
    !#################### PETSC DATA TYPES #########################
    PetscInt :: row, col(7),col1,ncolsmax
    PetscInt, Parameter :: i1=1,i2=2,i3=3,i5=5,i7=7
    PetscScalar :: val(7),val1,valq,valt,tol,minres
    Vec :: SOL_Vec,B_Vec
    Mat :: A_Mat
    !###############################################################

contains

subroutine distributeLoad(NLOCAL)

    use boundaryModule
    use indexModule
    use varModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    PetscErrorCode :: ierr

    integer, intent(in) :: NLOCAL
    Integer, Dimension(:), allocatable :: DNNZ(:),ONNZ(:)

    allocate(DNNZ(0:N-1),ONNZ(0:N-1))

    DNNZ=1
    ONNZ=0

    do B=1,NB
        call setBlockInd(B)

    !inner cells
        !print *, 'EAST'
        do K=2,NKM
        do I=2,NIM-1
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            IJKP=MIJK(IJK)
            !print *, IJKP
            DNNZ(IJKP)=DNNZ(IJKP)+1
            DNNZ(IJKP+NJCV)=DNNZ(IJKP+NJCV)+1
        end do
        end do
        end do
            
        !print *, 'NORTH'
        do K=2,NKM
        do I=2,NIM
        do J=2,NJM-1
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            IJKP=MIJK(IJK)
            !print *, IJKP
            DNNZ(IJKP)=DNNZ(IJKP)+1
            DNNZ(IJKP+1)=DNNZ(IJKP+1)+1
        end do
        end do
        end do

        !print *, 'TOP'
        do K=2,NKM-1
        do I=2,NIM
        do J=2,NJM
            IJK=IJKST+(K-1)*NI*NJ+(I-1)*NJ+J
            IJKP=MIJK(IJK)
            !print *, IJKP
            DNNZ(IJKP)=DNNZ(IJKP)+1
            DNNZ(IJKP+(NICV*NJCV))=DNNZ(IJKP+(NICV*NJCV))+1
        end do
        end do
        end do

    !block boundaries, identify local and offproc neighbours
        !print *, 'FACES: ',FACEST,NFACE
        do F=FACEST+1,FACEST+NFACE
            IJKP=MIJK(L(F))
            IJKB=MIJK(R(F))
            if (IJKB.ge.IJKPROC+N) then
                ONNZ(IJKP)=ONNZ(IJKP)+1
            else
                DNNZ(IJKP)=DNNZ(IJKP)+1
            end if
        end do
    end do


    ! Create Matrix
    call MatCreate(PETSC_COMM_WORLD,A_Mat,ierr)
    call MatSetSizes(A_Mat,N,N,PETSC_DECIDE,PETSC_DECIDE,ierr)
    call MatSetFromOptions(A_Mat,ierr)
    call MatMPIAIJSetPreallocation(A_Mat,PETSC_NULL_INTEGER,DNNZ,PETSC_NULL_INTEGER,ONNZ,ierr)
    call MatSeqAIJSetPreallocation(A_Mat,PETSC_NULL_INTEGER,DNNZ,ierr)
    call MatSetUp(A_Mat,ierr) !not necessary because already using Preallocation routine

    ! Create Vectors
    call VecCreate(PETSC_COMM_WORLD,SOL_Vec,ierr)
    call VecSetSizes(SOL_Vec,N,PETSC_DECIDE,ierr)
    call VecSetFromOptions(SOL_Vec,ierr)
    call VecDuplicate(SOL_Vec,B_Vec,ierr)

    call VecCreate(PETSC_COMM_WORLD,DTX_vec,ierr)
    call VecSetSizes(DTX_Vec,NXYZA,PETSC_DECIDE,ierr)
    call VecSetFromOptions(DTX_Vec,ierr)
    call VecDuplicate(DTX_Vec,DTY_Vec,ierr)
    call VecDuplicate(DTX_Vec,DTZ_Vec,ierr)

    call VecCreateSeq(PETSC_COMM_WORLD,NFACEAL,TR_Vec,ierr)
    call VecDuplicate(TR_Vec,DTXR_Vec,ierr)
    call VecDuplicate(TR_Vec,DTYR_Vec,ierr)
    call VecDuplicate(TR_Vec,DTYR_Vec,ierr)

end subroutine distributeLoad

end module coefModule

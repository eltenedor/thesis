module coefModule

    use parameterModule
    implicit none
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>


    real(KIND=PREC) :: Q(NXYZA), AP(NXYZA), AS(NXYZA), AN(NXYZA),&
                       AW(NXYZA), AE(NXYZA), AB(NXYZA), AT(NXYZA),&
                       AF(NFACEAL)
                       
    Mat :: Amat
    Vec :: solvec,bvec
    PetscInt, Parameter :: i1=1,i2=2,i3=3,i5=5,i7=7
    PetscInt :: row, col(7),col1,ncolsmax
    PetscScalar :: val(7),val1,valq,valt,tol,minres
    Integer, Dimension(:), allocatable :: DNNZ(:),ONNZ(:)

contains

subroutine distributeLoad(NLOCAL)

    use boundaryModule
    use indexModule
    implicit none
!#include <finclude/petscvec.h>
!#include <flinclude/petscmat.h>

    integer, intent(in) :: NLOCAL

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

    print *, 'DONE'

end subroutine distributeLoad

end module coefModule

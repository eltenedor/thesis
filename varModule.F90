! this module contains the arrays for storing the scalar T, its respective
! cell center gradient components and the mass fluxes F1...3
module varModule

    use parameterModule
    implicit none
#include <finclude/petscvecdef.h>

    real(KIND=PREC) ::  T(NXYZA),TO(NXYZA),TR(NFACEAL),TIME,&
                        DTX(NXYZA),DTY(NXYZA),DTZ(NXYZA),&
                        DTXR(NFACEAL),DTYR(NFACEAL),DTZR(NFACEAL),&
                        F1(NXYZA),F2(NXYZA),F3(NXYZA),&
                        FDIR(NDIRAL),FNEU(NNEUAL)

    Vec :: DTX_Vec,DTY_Vec,DTZ_Vec,DTXR_Vec,DTYR_Vec,DTZR_Vec,TR_Vec

contains

subroutine VecToArr(N,MAP,INP_Vec,OUTP)

    use parameterModule
    implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
    
    PetscErrorCode :: ierr
    Vec, intent(inout) :: INP_Vec
    Vec :: TEMP_Vec
    PetscScalar, pointer :: TEMP_Sca(:)
    PetscInt :: start_int,endvec_int
    VecScatter :: scatter
    IS :: FROM_Is, TO_Is

    integer, intent(in) :: N,MAP(N)
    integer :: I
    real(kind=PREC) ,intent(inout) :: OUTP(N)

    call VecCreateSeq(PETSC_COMM_SELF,N,TEMP_Vec,ierr)
    call ISCreateGeneral(PETSC_COMM_SELF,N,MAP,PETSC_COPY_VALUES,FROM_Is,ierr)
    call ISCreateStride(PETSC_COMM_SELF,N,0,1,TO_Is,ierr)
    call VecScatterCreate(INP_Vec,FROM_Is,TEMP_Vec,TO_Is,scatter,ierr)
    call VecScatterBegin(scatter,INP_Vec,TEMP_Vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecScatterEnd(scatter,INP_Vec,TEMP_Vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
    call VecGetArrayF90(TEMP_Vec,TEMP_Sca,ierr)

    do I=1,N
        OUTP(I)=TEMP_Sca(I)
    end do

    call VecRestoreArrayF90(TEMP_Vec,TEMP_Sca,ierr)

end subroutine VecToArr

end module varModule

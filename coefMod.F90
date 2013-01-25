module coef

    use param
    implicit none

#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    real(KIND=PREC) :: Q(NXYA), AP(NXYA), AS(NXYA), AN(NXYA),&
                       AW(NXYA), AE(NXYA)
                       
    Mat :: A
    Vec :: sol,b
    PetscInt, Parameter :: i1=1,i2=2,i5=5
    PetscInt :: row, col(5)
    PetscScalar :: val(5),valq,valt

end module coef

module coef

    use param
    implicit none
#include <finclude/petscvecdef.h>
#include <finclude/petscmatdef.h>


    real(KIND=PREC) :: Q(NXYA), AP(NXYA), AS(NXYA), AN(NXYA),&
                       AW(NXYA), AE(NXYA), AB(NXYZA), AT(NXYZA)
                       
    Mat :: A
    Vec :: sol,b
    PetscInt, Parameter :: i1=1,i2=2,i5=5,i7=7
    PetscInt :: row, col(5)
    PetscScalar :: val(5),valq,valt,tol,minres

end module coef

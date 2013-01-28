module coef

    use param
    implicit none

#include <finclude/petscvec.h>
#include <finclude/petscmat.h>

    real(KIND=PREC) :: Q(NXYA), AP(NXYA), AS(NXYA), AN(NXYA),&
                       AW(NXYA), AE(NXYA)
                       
    Mat :: A,A2
    Vec :: sol,b,vt1,vt2,res
    PetscInt, Parameter :: i1=1,i2=2,i5=5
    PetscInt :: row, col(5)
    PetscScalar :: val(5),valq,valt,tol,minres

end module coef
